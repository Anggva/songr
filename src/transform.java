import java.util.Arrays;
import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.DataLine;
import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.File;

public class transform {
        private int n, m;
        // Lookup tables. Only need to recompute when size of FFT changes.
        private double[] cos;
        private double[] sin;
        private double[] window;
        public transform(int n)
        {
                this.n = n;
                this.m = (int) (Math.log(n) / Math.log(2));
                // Make sure n is a power of 2
                if (n != (1 << m))
                        throw new RuntimeException("FFT length must be power of 2");
                // precompute tables
                cos = new double[n / 2];
                sin = new double[n / 2];
                for (int i = 0; i < n / 2; i++)
                {
                        cos[i] = Math.cos(-2 * Math.PI * i / n);
                        sin[i] = Math.sin(-2 * Math.PI * i / n);
                }
                makeWindow();
        }
    // JLS
    /** Makes a Hann (Hanning) window.
     * http://www.ece.uci.edu/docs/hspice/hspice_2001_2-220.html
     * http://blinkdagger.com/matlab/matlab-windowing-part-3/
     * http://en.wikipedia.org/wiki/Window_function
     */
    private void makeWindow()
        {
        window = new double[n];
        for (int i = 0; i < n; i++) {
            window[i] = 0.5 * (1 - Math.cos( (2 * Math.PI * i) / (n - 1) ));
        }
        }
    private double[] getWindow()
        {
                return window;
        }
    public double[] fftOfReals(float[] x) {
        double[] d = new double[x.length];
        for (int i = 0; i < d.length; i++) {
            d[i] = x[i];
        }
        return fftOfReals(d);
    }

    // JLS
    /** Given array x of n reals, returns an array of n/2 reals representing
     * the normalized magnitudes of the non-redundant complex values resulting from the
     * FFT of x. The magnitudes are normalized (divided) by n/2.
     * A windowing function is applied to x to prevent spectral leakage.
     * This function modifies x.
     * http://www.dsptutor.freeuk.com/analyser/SpectrumAnalyser.html (source code)
     */
    public double[] fftOfReals(double[] x) {
        // Apply windowing function
        double[] win = getWindow();
        for (int i = 0; i < x.length; i++) {
            x[i] = x[i] * win[i];
        }
        // Create array of imaginaries (all 0)
        double[] y = new double[x.length];
        Arrays.fill(y, 0.0);
        fft(x, y);
        // Take magnitudes of FFT
        double[] m = new double[x.length / 2];
        for (int i = 0; i < m.length; i++) {
            m[i] = Math.sqrt(x[i]*x[i] + y[i]*y[i]) / m.length;
        }
        return m;
    }

        /***************************************************************************
         * fft.c Douglas L. Jones University of Illinois at Urbana-Champaign January
         * 19, 1992 http://cnx.rice.edu/content/m12016/latest/
         *
         * fft: in-place radix-2 DIT DFT of a complex input
         *
         * input: n: length of FFT: must be a power of two m: n = 2**m input/output
         * x: double array of length n with real part of data y: double array of
         * length n with imag part of data
         *
         * Permission to copy and use this program is granted as long as this header
         * is included.
         **************************************************************************/
private void fft(double[] x, double[] y){
	int i, j, k, n1, n2, a;
    double c, s, t1, t2;
    // Bit-reverse
    j = 0;
    n2 = n / 2;
    for (i = 1; i < n - 1; i++){
    	n1 = n2;
        while (j >= n1){
        	j = j - n1;
            n1 = n1 / 2;
        }
        j = j + n1;
        if (i < j){
        	t1 = x[i];
        	x[i] = x[j];
        	x[j] = t1;
        	t1 = y[i];
        	y[i] = y[j];
        	y[j] = t1;
        }
    }
    //FFT
    n1 = 0;
    n2 = 1;
    for (i = 0; i < m; i++){
    	n1 = n2;
    	n2 = n2 + n2;
    	a = 0;
    	for (j = 0; j < n1; j++){
    		c = cos[a];
    		s = sin[a];
    		a += 1 << (m - i - 1);
    		for (k = j; k < n; k = k + n2){
    			t1 = c * x[k + n1] - s * y[k + n1];
    			t2 = s * x[k + n1] + c * y[k + n1];
    			x[k + n1] = x[k] - t1;
    			y[k + n1] = y[k] - t2;
    			x[k] = x[k] + t1;
    			y[k] = y[k] + t2;
    		}
    	}
    }
}
private int maxIndex (double [] magn){
	double max = magn[1];
	int k = 1;
	for (int i = 2; i < magn.length; i++)
		if (magn[i] > max){
			max = magn[i];
			k = i;
		}
	return k;
}
        // Test the FFT to make sure it's working
        public static void main(String[] args)
        {
                int N = 8;

                transform fft = new transform(N);
                File input = new File("Piano.pp.C4.wav");
                AudioInputStream audioStream;
                try{
                	audioStream = AudioSystem.getAudioInputStream(input);
                	byte[] bytes = new byte[(int) (audioStream.getFrameLength()) * (audioStream.getFormat().getFrameSize())]; 
                	audioStream.read(bytes);
                	int cSize = 4096;
                	int chunkN = bytes.length/cSize;
                	double avg = 0;
                	double [] magnFFT = new double [bytes.length];
                    for (int i = 0; i < chunkN; i++){
                    	double [] temp = new double [cSize];
                    	double[] im = new double[cSize];
                    	for (int j = 0; j < cSize; j++)
                    		temp[j] = bytes[cSize*i+j];
                    	beforeAfter(fft,temp,im);
                    	for (int j = 0; j < cSize; j++){
                    		magnFFT[cSize*i+j] = Math.sqrt(temp[j]*temp[j] + im[j]*im[j]);
                    		avg += magnFFT[cSize*i+j];
                    	} 
                    	avg /= cSize;
                    }
                    //Average magnitude in the array
                    avg /= chunkN;
                    
                    double [] nmfft = new double [chunkN * cSize];
                    for (int i = 0; i < nmfft.length; i++){
                    	nmfft[i] = magnFFT[i];
                    	//System.out.print(nmfft[i]+" ");
                    	//if (i%50==0)
                    		//System.out.println("");
                    }
                    int largN = 10;
                    //Searching for the largN largest magnitudes
                    double [] topMag = new double [largN];
                    int size = nmfft.length;
                    double topAvg= 0;
                    for (int i = 0; i < largN; i++){
                    	int j = fft.maxIndex(nmfft);
                    	topMag[i] = nmfft[j];
                    	topAvg += topMag[i];
                    	for (int k = j; k < size-1; k++)
                    		nmfft[k] = nmfft[k+1];
                    }
                    topAvg /= largN;
                    System.out.println(topAvg);
                    System.out.println(avg);

                }
                catch (Exception e){
                 e.printStackTrace();
                }
        }

        private static void beforeAfter(transform fft, double[] re, double[] im)
        {
                fft.fft(re, im);
        }

        private static void printReIm(double[] re, double[] im)
        {
                System.out.print("Re: [");
                for (int i = 0; i < re.length; i++)
                        System.out.print(((int) (re[i] * 1000) / 1000.0) + " ");

                System.out.print("]\nIm: [");
                for (int i = 0; i < im.length; i++)
                        System.out.print(((int) (im[i] * 1000) / 1000.0) + " ");

                System.out.println("]");
        }
}