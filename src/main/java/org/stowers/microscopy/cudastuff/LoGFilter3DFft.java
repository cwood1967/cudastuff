package org.stowers.microscopy.cudastuff;

import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import net.imagej.ImageJ;
import ij.ImagePlus;

import org.scijava.command.Command;
import org.scijava.command.Previewable;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.stowers.microscopy.utils.ExpLookup;

import edu.emory.mathcs.jtransforms.fft.FloatFFT_3D;

@Plugin(type = Command.class, name = "Click for Roi",  menuPath="Plugins>Chris>NotReady>LoG")
public class LoGFilter3DFft implements Previewable, Command {

	@Parameter
	ImagePlus imp;
	
	
	int sizeX;
	int sizeY;
	int sizeZ;
	
	@Parameter(label="Sigma X:")
	float sigX;
	
	@Parameter(label="Sigma Y:")
	float sigY;
	
	@Parameter(label="Sigma Z:")
	float sigZ;
	
	float sigX2, sigY2, sigZ2;
	
	float[] filter;
	float[] imagePixels;
	float[] filteredPixels;
	ImageStack stack;
	ImagePlus fimp;
	ImageStack fstack;
	ExpLookup exp;
	
	public LoGFilter3DFft() {
		
	}
	
	public LoGFilter3DFft(int sizeX, int sizeY, int sizeZ,
			float sigX, float sigY, float sigZ) {
	
		this.sizeX = sizeX;
		this.sizeY = sizeY;
		this.sizeZ = sizeZ;
		this.sigX = sigY;
		this.sigY = sigY;
		this.sigZ = sigZ;
		
		sigX2 = sigX*sigX;
		sigY2 = sigY*sigY;
		sigZ2 = sigZ*sigZ;
	
		exp = new ExpLookup();
		stack = new ImageStack(sizeX, sizeY, sizeZ);
	
		imp = new ImagePlus("LogStack");
		
	}
	
	@Override
	public void run() {
		
		sigX2 = sigX*sigX;
		sigY2 = sigY*sigY;
		sigZ2 = sigZ*sigZ;
		
		sizeX = imp.getWidth();
		sizeY = imp.getHeight();
		sizeZ = imp.getStackSize();
		
		exp = new ExpLookup();
		long t1 = System.currentTimeMillis();
		makePixelsArray();
		long t2 = System.currentTimeMillis();
		System.out.println(t2 -t1);
		calcLoG();
		filteredLoG();
		fimp = new ImagePlus("Filter Stack", fstack);
		fimp.show();
	}

	private float[] arrayMax(float[] a) {
		
		float max = -Float.MAX_VALUE;
		float min = Float.MAX_VALUE;
		
		for (int i = 0; i < a.length; i++) {
			if (a[i] > max) max = a[i];
			if (a[i] < min) min = a[i];
		}
		
		return new float[] {min, max};
	}
	public void filteredLoG() {
	
		FloatFFT_3D fft3d = new FloatFFT_3D(sizeZ, sizeY, sizeX);
		
		int np = 2*sizeX*sizeY*sizeZ;
		
		float[] ifilter = new float[np];
		
		System.arraycopy(filter, 0, ifilter, 0, np);
		
		float[] b = arrayMax(filter);
		System.out.println("Filter " + b[0] + " " + b[1] + " " + filter.length);
		b = arrayMax(imagePixels);
		System.out.println("imagePixels " + b[0] + " " + b[1] + " " + filter.length);
		
		b = arrayMax(ifilter);
		System.out.println("iFilter " + b[0] + " " + b[1] + " " + ifilter.length);
		
		fft3d.complexForward(ifilter);
		
		fft3d.complexForward(imagePixels);
		
		b = arrayMax(filter);
		System.out.println("Filter " + b[0] + " " + b[1] + " " + filter.length);
		
		b = arrayMax(ifilter);
		System.out.println("iFilter " + b[0] + " " + b[1] + " " + ifilter.length);
		
		for (int i = 0; i < ifilter.length; i++) {
			imagePixels[i] = imagePixels[i]*ifilter[i];
		}
		
		b = arrayMax(imagePixels);
		System.out.println("imagePixels " + b[0] + " " + b[1] + " " + ifilter.length);
		fft3d.complexInverse(imagePixels, false);
		
		b = arrayMax(imagePixels);
		System.out.println("imagePixels " + b[0] + " " + b[1] + " " + ifilter.length);
		
		filteredPixels = new float[imagePixels.length/2];
		int counter = 0;
		
		for (int i = 0; i < filteredPixels.length; i++) {
//			filteredPixels[i] = imagePixels[counter]*imagePixels[counter] +
//					imagePixels[counter + 1]*imagePixels[counter + 1];
			filteredPixels[i] = imagePixels[counter];
			counter += 2;
		}
		

//		float[] fpminmax = arrayMax(filteredPixels);
//		for (int i = 0; i < filteredPixels.length; i++) {
//			filteredPixels[i] = 200.f*filteredPixels[i]/(fpminmax[1] - fpminmax[0]);
//		}
		System.out.println("filteredPixels " + b[0] + " " + b[1] + " " + filteredPixels.length);
		fstack = new ImageStack(sizeX, sizeY, sizeZ);
		
		float[] shifted = new float[filteredPixels.length];
		
		int index = 0;
		byte quads = 15;
		int xoff = sizeX/2 - 1;
		int yoff = sizeY/2 - 1;
		int y;
		int x;
		for (int k = 0; k < sizeZ; k++) {
			float[] px = new float[sizeX*sizeY];
			int ns = 0;
			int size = sizeX*sizeY;
			
			for (int j = 0; j < sizeY; j++) {
			
				if (j < sizeY/2) {
					y = j + yoff;
				} else {
					y = j - yoff;
				}
				
				for (int i = 0; i < sizeX; i++) {
					if (i < sizeX/2) {
						x = i + xoff;
					} else {
						x = i - xoff;
					}
					
					int pIndex = sizeX*y + x;
					try {
						px[pIndex] = filteredPixels[index];
					}
					catch (Exception e) {
						System.out.println(x + " " + y + " " + pIndex + " " + index) ;
						e.printStackTrace();
					}
					index++;
				}
				
			}
			fstack.setPixels(px,  k + 1);
		}
	}
	
	private void makePixelsArray() {
		imagePixels = new float[2*sizeX*sizeY*sizeZ];
		
		ImageStack pstack = imp.getImageStack();
		float[] p;
		
		int index = 0;
		for (int k = 0; k < sizeZ; k++) {
			p = (float[])pstack.getProcessor(k + 1).convertToFloatProcessor().getPixels();
			for (int i = 0; i < p.length; i++) {
				imagePixels[index] = p[i];
				index += 2;
			}
		}
	}
	
	protected void calcLoG() {
		
		float x0 = sizeX/2.f;
		float y0 = sizeY/2.f;
		float z0 = sizeZ/2.f;
		float x, y, z;
		
		long t1 = System.currentTimeMillis();
		float xcof, ycof, zcof;
		float xkern, ykern, zkern;
		float g;
		
		float[][] pixels = new float[sizeZ][sizeX*sizeY];
		int index;
		
		float isigZ2 = 1.f/sigZ2;
		float isigY2 = 1.f/sigY2;
		float isigX2 = 1.f/sigX2;
		
		float minKern = 0;
		filter = new float[sizeX*sizeY*2*sizeZ];  //notice the *2 in the middle
		float maxf = 0;
		ImageStack xstack = new ImageStack(sizeX, sizeY, sizeZ);
		for (int k =0; k < sizeZ; k++) {
			z = (k - z0);
			zcof = (z*z - sigZ2)/(isigZ2*isigZ2);
			zkern = 0.5f*z*z*(isigZ2);
			
			index = 0;
			
			for (int j = 0; j < sizeY; j++) {
				y = j - y0;
				ycof = (y*y - sigY2)/(isigY2*isigY2);
				ykern = 0.5f*y*y*(isigY2);	
				
				for (int i = 0; i < sizeX; i++) {
					x = i - x0;
					xcof = (x*x - sigX2)/(isigX2*isigX2);
					xkern = 0.5f*x*x*(isigX2);
					
//					g = (float)Math.exp(-(xkern + ykern + zkern));
					g = (float)exp.texp(-(xkern + ykern + zkern));
					g *= -(xcof + ycof + zcof);
//					
					pixels[k][index/2] = g;
					
					if (g > maxf) {
						maxf = g;
					}
					filter[index] = g;
					index += 2;
				}
				
				xstack.setPixels(pixels[k],  k + 1);
			}
			
		}
		
//		float[] filter = new float[sizeX*sizeY*sizeZ];
		long t2 = System.currentTimeMillis();
		System.out.println("Time: " + .001*(t2 -t1));
		System.out.println("Max f: " + maxf);
		ImagePlus ximp = new ImagePlus("Filter");
		ximp.setStack(xstack);
		ximp.show();
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

		final ImageJ imagej = net.imagej.Main.launch(args);
		
//		LoGFilter3DFft log = new LoGFilter3DFft(1024, 1024, 64, 8, 8, 4);
//		
//		log.calcLoG();
//		
//		long t1 = System.currentTimeMillis();
//		log.filteredLoG();
//		long t2 = System.currentTimeMillis();
//		System.out.println(t2 - t1);
		
	}

	
	@Override
	public void cancel() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void preview() {
		// TODO Auto-generated method stub
		
	}

}
