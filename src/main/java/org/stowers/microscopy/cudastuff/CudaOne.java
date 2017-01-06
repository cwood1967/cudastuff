package org.stowers.microscopy.cudastuff;

import ij.IJ;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;

import static jcuda.jcufft.JCufft.CUFFT_FORWARD;
import static jcuda.jcufft.JCufft.cufftDestroy;
import static jcuda.jcufft.JCufft.cufftExecC2C;
import jcuda.jcufft.JCufft;

import static jcuda.jcufft.JCufft.cufftPlan2d;

import org.scijava.command.Command;
import org.scijava.command.Previewable;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import jcuda.jcufft.cufftHandle;
import jcuda.jcufft.cufftType;
import net.imagej.patcher.LegacyInjector;
import net.imagej.ImageJ;

import edu.emory.mathcs.jtransforms.fft.FloatFFT_2D;

import java.util.Random;

@Plugin(type = Command.class, name = "Click for Roi",  menuPath="Plugins>Chris>Cufft")
public class CudaOne implements  Previewable, Command {

	static {
        LegacyInjector.preinit();
    }
	
	@Parameter
	ImagePlus imp;
	
	FloatProcessor fip;
	
	public void run() {
		
		fip = imp.getProcessor().convertToFloatProcessor();
//		fip = new FloatProcessor(128,  128);
//		fip.setPixels(makesine(128));
//		ImagePlus imp = new ImagePlus("Input", fip);
//		imp.show();
//		IJ.log("Working on FFT");
		doFFT();
	}
	
	private void doFFT() {
		
		int w = fip.getWidth();
		int h = fip.getHeight();
		
		float[] tmp = (float[])fip.getPixels();
		
		int[] size = padSize(w, h);
		float[] padded = pad(tmp, w, h, size[0], size[1]);
		float[] pixels = new float[2*padded.length];
		
		Random rand = new Random();
		for (int i = 0; i < padded.length; i += 2) {
			pixels[i] = padded[i];
			pixels[i + 1] = 0; //(float)(2.*rand.nextGaussian());
		}
 		
		float[] fixels = new float[pixels.length]; //[h*(w/2 +1) - 20];
		
		float[] jpixels = new float[padded.length];
		System.arraycopy(padded,0, jpixels, 0, padded.length); 
		FloatFFT_2D jfft = new FloatFFT_2D(size[1], size[0]);
		jfft.realForward(jpixels);
		
		ImageProcessor jfftip = new FloatProcessor(size[0],  size[1]);
		float[] jfftPixels = new float[jpixels.length/2];
		jfftip.setPixels(jfftPixels);
		ImagePlus jfftimp = new ImagePlus("jtransforms fft");
		jfftimp.setProcessor(jfftip);
		
		
		cufftHandle plan = new cufftHandle();
		cufftPlan2d(plan, size[1], size[0], cufftType.CUFFT_C2C);
		JCufft.cufftExecC2C(plan, pixels, fixels, JCufft.CUFFT_FORWARD);
		JCufft.cufftDestroy(plan);
		float[] c = new float[pixels.length/2];
		float[] ps = new float[pixels.length/2];
		float[] p = new float[pixels.length/2];
		int j = 0;
		for (int i = 0; i < p.length; i++) {
			if (i < 20) {
				System.out.println(i + " " + j + " " + fixels[j] + " " + fixels[j + 1]);
			}
			
			ps[i] = (float)(fixels[j]*fixels[j] + fixels[j + 1]*fixels[j + 1]);
			if (i < 65536 && j < 65535) {
				jfftPixels[i] = (float)(jpixels[j]*jpixels[j] + jpixels[j + 1]*jpixels[j + 1]);
			}
			p[i] = jpixels[j];
		
			c[i] = jpixels[j + 1];
			j += 2;
		}
//		c = fixels;
		System.out.println(j + " " + fixels.length);
		int nw = size[0];
		int nh = size[1];
		ImageProcessor xip = new FloatProcessor(nw,  nh);
		xip.setPixels(p);
		ImagePlus xmp = new ImagePlus("FFT Real", xip);
		xmp.show();
		 
		ImageProcessor cip = new FloatProcessor(nw,  nh);
		cip.setPixels(c);
		ImagePlus cmp = new ImagePlus("FFT Imag", cip);
		cmp.show();
		
		ImageProcessor psip = new FloatProcessor(nw,  nh);
		psip.setPixels(ps);
		ImagePlus psmp = new ImagePlus("FFT ps", psip);
		psmp.show();
		
		ImageProcessor padip = new FloatProcessor(nw,  nh);
		padip.setPixels(padded);
		ImagePlus padmp = new ImagePlus("Padded", padip);
		padmp.show();
		
		jfftimp.show();
	}
	
	private int[] padSize(int w, int h) {
		int nw = (int)(Math.ceil(Math.log(w)/Math.log(2.)));
		int nh = (int)(Math.ceil(Math.log(h)/Math.log(2.)));
		
		nw = (int)Math.pow(2, nw);
		nh = (int)Math.pow(2, nh);
		
		if (nw >= nh) {
			nh = nw;
		} else {
			nw =nh;
		}
		return new int[] {nw, nh};
	}
	private float[] pad(float[] pixels, int w, int h, int nw, int nh) {
		
		float[] res = new float[nw*nh];
		float v = (float)fip.getStatistics().mean;
		for (int i =0; i < res.length; i++) {
			res[i] = v;
		}
		
		int woffset = 0*(nw - w)/2;
		int hoffset = 0*(nh - h)/2;
		for (int j = 0; j < h; j++) {
			int y = hoffset + j;
			for (int i = 0; i < w; i++) {
				int x = woffset + i;
				int index = y*nw + x;
				int pindex = j*w + i;
				res[index] = pixels[pindex];
				
			}
		}
		
		return res;
	}
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
	        final ImageJ imagej = net.imagej.Main.launch(args);
	}

	@Override
	public void cancel() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void preview() {
		// TODO Auto-generated method stub
		
	}
	
	public float[] makelines(int size) {
		float[] res = new float[size*size];
		
		int ix = size/2;
		for (int iy = 0; iy < size; iy++) {
			int index = iy*size + ix;
			res[index] = 0;
		}
		
		int iy = size/2;
		for (ix = 0; ix < size; ix++) {
			int index = iy*size + ix;
			res[index] = 20;
		}
		
		return res;
	}
	public float[] makefunc(int size) {
		
		float[] res = new float[size*size];
		
		double ystart = -6.;
		double ystop = 6;
		double xstart = -6;
		double xstop = 6;
		double del = (xstop - xstart)/size;
		 
		for (int iy = 0; iy < size; iy++) {
			double y = ystart + iy*del;
			double fy = -y*y/2.;
			for (int ix = 0; ix < size; ix++) {
				int index = iy*size + ix;
				double x = xstart + ix*del;
				double fx = -x*x/2.;
				double f1 = Math.exp(fx + fy);
				res[index] = (float)f1;
			}
		}
		return res;
	}
	
public float[] makesine(int size) {
		
		float[] res = new float[size*size];
		
		double ystart = 0.;
		double ystop = 2.*Math.PI;
		double xstart = 0;
		double xstop = 2.*Math.PI;
		double del = (xstop - xstart)/size;
		 
		for (int iy = 0; iy < size; iy++) {
			double y = ystart + iy*del;
			double fy = -y*y/2.;
			for (int ix = 0; ix < size; ix++) {
				int index = iy*size + ix;
				double x = xstart + ix*del;
				double fx = -x*x/2.;
				double f1 = Math.sin(2.*x);
				res[index] = (float)f1;
			}
		}
		return res;
	}

}
