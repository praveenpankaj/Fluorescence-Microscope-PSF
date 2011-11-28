package plugins.praveen.PSF;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_2D;
//import edu.emory.mathcs.jtransforms.fft.FloatFFT_2D;
//import loci.poi.util.ByteField;
import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import icy.type.DataType;
//import icy.type.TypeUtil;

//import plugins.adufour.ezplug.EzVarDouble;
//import plugins.adufour.ezplug.EzVarInteger;

public class PSFCalculator {
	private int _w;
	private int _h;
	private int _z;
	private double _indexRefr;
	private double _na;
	private int _lex;
	private int _lem;
	private double _sa;
	private double _xySampling;
	private double _zSampling;
	private double _depth;
	
	public final static int DEFAULT_W = 512;
	public final static int DEFAULT_H = 512;
	public final static int DEFAULT_Z = 64;
	public final static double DEFAULT_INDEXREFR = 1.515;
	public final static double DEFAULT_NA = 0.5;
	public final static int DEFAULT_LEX = 488;
	public final static int DEFAULT_LEM = 520;
	public final static double DEFAULT_INDEXSP = 1.44;
	public final static double DEFAULT_XYSAMPLING = 50.00;
	public final static double DEFAULT_ZSAMPLING = 200.00;
	public final static double DEFAULT_DEPTH = 0.0;
	PSFCalculator() {
		setW(DEFAULT_W);
		setH(DEFAULT_H);
		setZ(DEFAULT_Z);
		setXYSAMPLING(DEFAULT_XYSAMPLING);
		setZSAMPLING(DEFAULT_ZSAMPLING);
		setIndexRefr(DEFAULT_INDEXREFR);
		setNA(DEFAULT_NA);
		setLEX(DEFAULT_LEX);
		setLEM(DEFAULT_LEM);
		setSA(DEFAULT_INDEXSP);
		
		setDEPTH(DEFAULT_DEPTH);		
	}
	public Sequence compute(){
		
		final DoubleFFT_2D fft = new DoubleFFT_2D(_w, _h);
		//IcyBufferedImage psf3d = new IcyBufferedImage(_w, _h, 1, DataType.DOUBLE);
		
		int hc = _h/2;
		int wc = _w/2;
		int zc = _z/2;
		double kSampling = (2*Math.PI)/(_h*_xySampling); //Fourier space sampling
		double lambdaObj = _lem/_indexRefr;//Wavelength of light inside the medium
		double k0 = (2*Math.PI)/_lem;//Wave vector
		double kObj = (2*Math.PI)/lambdaObj;//Wave vector in the medium
		double kMax = (2*Math.PI*_na)/(_lem*kSampling);//Maximum aperture radius
		double saCoeff = _depth * _sa;
		
		// Define the zero defocus pupil function
		IcyBufferedImage pupil = new IcyBufferedImage(_w, _h, 2, DataType.FLOAT); // channel 1 is real and channel 2 is imaginary
		float[] pupilRealBuffer = pupil.getDataXYAsFloat(0);//Real
		float[] pupilImagBuffer = pupil.getDataXYAsFloat(1);//imaginary
		
		IcyBufferedImage dpupil = new IcyBufferedImage(_w, _h, 2, DataType.DOUBLE); // channel 1 is real and channel 2 is imaginary
		double[] dpupilRealBuffer = dpupil.getDataXYAsDouble(0); //Real
		double[] dpupilImagBuffer = dpupil.getDataXYAsDouble(1); //imaginary
		
		//Calculate the cosine and the sine components
		IcyBufferedImage ctheta = new IcyBufferedImage(_w, _h, 1, DataType.DOUBLE);
		double[] cthetaBuffer = ctheta.getDataXYAsDouble(0);
		IcyBufferedImage stheta = new IcyBufferedImage(_w, _h, 1, DataType.DOUBLE);
		double[] sthetaBuffer = stheta.getDataXYAsDouble(0);
		
		Sequence psf3d = new Sequence();
		psf3d.setName("double image");
        
		for (int k =  0 ; k < _z; k++)
    	{// Define the defocus pupils			
			
			double defocus = k-zc;
			defocus = defocus*_zSampling;	
			
			for(int x = 0; x < _w; x++)
			{
				for(int y = 0; y < _h; y++)
				{   
					double kxy = Math.sqrt( Math.pow(x-wc, 2) + Math.pow(y-hc, 2) );
        		
					pupilRealBuffer[pupil.getOffset(x, y)] = ((kxy < kMax) ? 1 : 0); //Pupil bandwidth constraints
					pupilImagBuffer[pupil.getOffset(x, y)] = 0; //Zero phase 
        		
					sthetaBuffer[x + y * _h] = Math.sin( kxy * kSampling / kObj );
					sthetaBuffer[x + y * _h] = (sthetaBuffer[x + y * _h]< 0) ? 0: sthetaBuffer[x + y * _h];
					cthetaBuffer[x + y * _h] = Double.MIN_VALUE + Math.sqrt(1 - Math.pow(sthetaBuffer[x + y * _h], 2));
					dpupilRealBuffer[x + y * _h] = pupilRealBuffer[pupil.getOffset(x, y)] * Math.cos((defocus * kObj*cthetaBuffer[x + y * _h])+saCoeff);
					dpupilImagBuffer[x + y * _h] = pupilRealBuffer[pupil.getOffset(x, y)] * Math.sin((defocus * kObj*cthetaBuffer[x + y * _h])+saCoeff);
        		}
			}
			double[] psf2d = dpupil.getDataCopyCXYAsDouble();
			fft.complexInverse(psf2d, false);
			
			IcyBufferedImage timg = new IcyBufferedImage(_w, _h, 1, DataType.DOUBLE);
			timg.beginUpdate();
			try{
				for(int x = 0; x < (wc+1); x++)
				{
					for(int y = 0; y < (hc+1); y++)
					{
						timg.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(psf2d[(((wc-x) + (hc-y) * _h)*2)+0],2)+Math.pow(psf2d[(((wc-x) + (hc-y) * _h)*2)+1], 2)));
						//timg.setDataAsDouble(x, y, 1, psf2d[(((wc-x) + (hc-y) * _h)*2)+1]);

					}
					for(int y = hc+1; y < _h; y++)
					{
						timg.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(psf2d[(((wc-x) + (_h+hc-y) * _h)*2)+0], 2)+Math.pow(psf2d[(((wc-x) + (_h+hc-y) * _h)*2)+1], 2)));
						//timg.setDataAsDouble(x, y, 1, psf2d[(((wc-x) + (_h+hc-y) * _h)*2)+1]);
					}
					
				}
				for(int x = (wc+1); x < _w; x++)
				{
					for(int y = 0; y < (hc+1); y++)
					{
						timg.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(psf2d[(((_w+wc-x) + (hc-y) * _h)*2)+0], 2)+Math.pow(psf2d[(((_w+wc-x) + (hc-y) * _h)*2)+1], 2)));
						//timg.setDataAsDouble(x, y, 1, psf2d[(((_w+wc-x) + (hc-y) * _h)*2)+1]);
					}
					for(int y = hc+1; y < _h; y++)
					{
						timg.setDataAsDouble(x, y, 0, Math.sqrt(Math.pow(psf2d[(((_w+wc-x) + (_h+hc-y) * _h)*2)+0], 2)+Math.pow(psf2d[(((_w+wc-x) + (_h+hc-y) * _h)*2)+1],2)));
						//timg.setDataAsDouble(x, y, 1, psf2d[(((_w+wc-x) + (_h+hc-y) * _h)*2)+1]);
					}
				}

			}finally {
			timg.endUpdate();
			}
			
            psf3d.addImage(timg);
    	}
		
	return psf3d;

	}
	public int getW() {
		return _w;
	}
	public int getH() {
		return _h;
	}
	public int getZ() {
		return _z;
	}
	public double getIndexRefr() {
		return _indexRefr;
	}
	public double getNA() {
		return _na;
	}
	public int getLEX() {
		return _lex;
	}
	public int getLEM() {
		return _lem;
	}
	public double getSA() {
		return _sa;
	}
	public double getXYSAMPLING() {
		return _xySampling;
	}
	public double getZSAMPLING() {
		return _zSampling;
	}
	public double getDEPTH() {
		return _depth;
	}
	
	public void setW(int src) {
		_w = src;
	}
	public void setH(int src) {
		_h = src;
	}
	public void setZ(int src) {
		_z = src;
	}
	public void setIndexRefr(double src) {
		_indexRefr = src;
	}
	public void setNA(double src) {
		_na = src;
	}
	public void setLEX(int src) {
		_lex = src;
	}
	public void setLEM(int src) {
		_lem = src;
	}
	public void setSA(double src) {
		_sa = src;
	}
	public void setXYSAMPLING(double src) {
		_xySampling = src;
	}
	public void setZSAMPLING(double src) {
		_zSampling = src;
	}
	public void setDEPTH(double src) {
		_depth = src;
	}
}
