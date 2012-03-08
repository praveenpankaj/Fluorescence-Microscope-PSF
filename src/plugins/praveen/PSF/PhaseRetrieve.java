package plugins.praveen.PSF;

import icy.gui.frame.progress.AnnounceFrame;
import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import icy.type.DataType;
import icy.type.collection.array.Array1DUtil;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.filtering.Convolution1D;
import plugins.adufour.filtering.Kernels1D;
import icy.math.ArrayMath;
import icy.math.MathUtil;

import javax.media.jai.BorderExtender;
//import javax.media.jai.JAI;
//import javax.media.jai.KernelJAI;
import javax.media.jai.RenderedOp;
import javax.media.jai.operator.BorderDescriptor;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_2D;


public class PhaseRetrieve extends EzPlug {
	EzVarSequence _input = new EzVarSequence("Choose the 3D fluorescence bead image");
	EzVarDouble _xySampling = new EzVarDouble("Image pixel spacing, in nm", 92.00, 10.00, 50000.00, 0.1);
	EzVarDouble _zSampling = new EzVarDouble("Slice spacing (z), in nm", 277.00, 10.00, 50000.00, 0.1);
	EzVarDouble _objNA = new EzVarDouble("Effective numerical aperture of the objective lens", 1.4, 0.001, 4.00, 0.01);
	EzVarDouble _indexImmersion = new EzVarDouble("Refractive index of the lens immersion medium", 1.518, 1.00, 4.00, 0.01);
	EzVarInteger _lem = new EzVarInteger("Emission peak wavelength, in nm", 520, 405, 750, 1);
	EzVarInteger _bgd = new EzVarInteger("Mean background fluorescence intensity", 0, 0, 100000, 1);
	EzVarDouble _sigma = new EzVarDouble("Gaussian filter parameter", 0.5, 0.1, 5.00, 0.01);
	EzVarDouble _alpha = new EzVarDouble("Step size for the iterative algorithm", 0.6, 0.5, 0.99, 0.01); 
	EzVarInteger _nIter = new EzVarInteger("Number of iterations", 30, 3, 10000, 1);


	@Override
	protected void initialize() {		

		super.addEzComponent(_input);
		super.addEzComponent(_xySampling);
		super.addEzComponent(_zSampling);
		super.addEzComponent(_objNA);
		super.addEzComponent(_indexImmersion);
		super.addEzComponent(_lem);
		super.addEzComponent(_bgd);
		super.addEzComponent(_sigma);
		super.addEzComponent(_alpha);
		super.addEzComponent(_nIter);		         
	}

	@Override
	protected void execute() {
		Sequence pupil = null;
		pupil = estimatepupil(_input.getValue(), _xySampling.getValue(), _zSampling.getValue(), _objNA.getValue(), _indexImmersion.getValue(), _lem.getValue(), _bgd.getValue(), _sigma.getValue(), _alpha.getValue(), _nIter.getValue());
		/*addSequence(pupil);		
		pupil.setName("Estimated Back Aperture Pupil");
		pupil.setChannelName(0, "Magnitude");
		pupil.setChannelName(1, "Phase");*/
		//MessageDialog.showDialog("Test is working fine!");
	}

	public Sequence estimatepupil(Sequence sequence, double _xySampling, double _zSampling, double _objNA, double _indexImmersion, int _lem, int _bgd, double _sigma, double _alpha, int _nIter) 
	{
		// TODO Auto-generated method stub
		Sequence pupil = new Sequence();
		pupil.setName("Estimated Pupil");
		pupil.setChannelName(0, "Amplitude");
		pupil.setChannelName(1, "Phase");
		Sequence psf3d = new Sequence();
		psf3d.setName("Estimated PSF");
		Sequence resizedSeq = new Sequence();
		final int NSECTIONS = 4;
		final int NRINGS = 3;

		int _w = sequence.getSizeX();
		int _h = sequence.getSizeY();
		int _z = sequence.getSizeZ();		

		//1. Calculate the parameter necessary for the algorithm
		double _lambdaObj = _lem/_indexImmersion; //Wavelength inside the lens
		double _kObj = (2 * Math.PI)/_lambdaObj; //Wavenumber inside the lens
		double _rMax = (NRINGS * 0.061 * _lambdaObj)/(_objNA * _xySampling);//Maximum spread
		double _k0 = (2*Math.PI)/_lem;//Wave vector

		// 1. Resize the input data to make it a square image
		int leftPad = 0;
		int rightPad = 0;
		int topPad = 0;
		int botPad = 0;
		if(_w>_h)
		{ 
			double dh = _w-_h;
			_h = _w;
			if(Math.IEEEremainder(dh, 2)==0)
			{
				topPad = (int)dh/2;
				botPad = (int)dh/2;

			}
			else
			{
				topPad = (int) Math.ceil(dh/2);
				botPad = (int) Math.floor(dh/2);
			}
		}
		else
		{
			double dw = _h -_w;
			_w = _h;
			if(Math.IEEEremainder(dw, 2)==0)
			{
				leftPad = (int)dw/2;
				rightPad = (int)dw/2;

			}
			else
			{
				leftPad = (int) Math.ceil(dw/2);
				rightPad = (int) Math.floor(dw/2);
			}
		}
		for(int iz=0;iz<_z;iz++)
		{
			//origImage = sequence.getImage(0, iz, 0).getScaledCopy(_w, _w, false, SwingConstants.CENTER, SwingConstants.CENTER);
			final RenderedOp renderedOp = BorderDescriptor.create(sequence.getImage(0, iz, 0), leftPad, rightPad, topPad, botPad, BorderExtender.createInstance(BorderExtender.BORDER_REFLECT), null);
			IcyBufferedImage resizedImage = IcyBufferedImage.createFrom(renderedOp.getAsBufferedImage());
			resizedImage.dataChanged();
			resizedSeq.addImage(resizedImage);
		}
		final DoubleFFT_2D fftOp = new DoubleFFT_2D(_w, _h);
		double kSampling = 2*Math.PI/(_w*_xySampling);
		double _kMax = (2 * Math.PI * _objNA)/(_lambdaObj*kSampling); //Maximum permissible frequency
		int hc = (int) Math.ceil(_h/2);
		int wc = (int) Math.ceil(_w/2);
		double[][] seqArray = new double[_z][_w*_h];
		double[][] bgRemovedArray = new double[_z][_w*_h];

		//2. Find central plane
		//Sequence zMaxProj = new Sequence();
		int cPlane = 0;
		double[] zMaxIntensity = new double[_z];
		double maxIntensity=0;
		//Sequence zMaxProj = Projection.zProjection(, Projection.ProjectionType.MAX, true);

		// 3. Remove Mean Background FLuorescence Intensity
		for(int iz = 0;iz<_z;iz++)
		{		
			IcyBufferedImage zImage = sequence.getImage(0, iz, 0);			
			zImage.updateComponentsBounds(true, true);
			zMaxIntensity[iz] = zImage.getComponentUserMaxValue(0);
			if(maxIntensity < zMaxIntensity[iz])
			{
				cPlane = iz;
				maxIntensity = zMaxIntensity[iz];
			}
			seqArray[iz] = Array1DUtil.arrayToDoubleArray(resizedSeq.getDataXY(0, iz, 0), false);	
			for(int ix = 0;ix<_w;ix++)
			{
				for(int iy = 0;iy<_h;iy++)
				{
					bgRemovedArray[iz][ix + iy*_w] = seqArray[iz][ix + iy*_w]-_bgd;
					bgRemovedArray[iz][ix + iy*_w] = ((bgRemovedArray[iz][ix + iy*_w] < 0) ? 0 : bgRemovedArray[iz][ix + iy*_w]);					
				}
			}			
		}
		cPlane = cPlane+1;
		//ArrayMath.divide(bgRemovedArray, ArrayMath.max(Array1DUtil.arrayToDoubleArray(bgRemovedArray, true)));

		//4. Display the focal plane information
		new AnnounceFrame("Detected focal plane at the " + cPlane + "th slice.");
		int[] selectedPlanes = new int[]{cPlane-15, cPlane-2, cPlane+2, cPlane+15};
		double[] defocus = new double[NSECTIONS];
		ArrayMath.subtract(selectedPlanes, cPlane, selectedPlanes);
		ArrayMath.multiply(Array1DUtil.arrayToDoubleArray(selectedPlanes, true), (double)_zSampling, defocus);	

		//5. Initialize Pupil Function
		// Define the zero defocus pupil function
		IcyBufferedImage pupilImage = new IcyBufferedImage(_w, _h, 2, DataType.DOUBLE); // channel 1 is magnitude and channel 2 is phase
		double[] pupilMagBuffer = pupilImage.getDataXYAsDouble(0);//Real
		double[] pupilPhBuffer = pupilImage.getDataXYAsDouble(1);//imaginary

		IcyBufferedImage dpupilImage = new IcyBufferedImage(_w, _h, 2, DataType.DOUBLE); // channel 1 is real and channel 2 is imaginary
		double[] dpupilReBuffer = dpupilImage.getDataXYAsDouble(0); //Real
		double[] dpupilImBuffer = dpupilImage.getDataXYAsDouble(1); //imaginary

		//6. Calculate the cosine and the sine components
		IcyBufferedImage ctheta = new IcyBufferedImage(_w, _h, 1, DataType.DOUBLE);
		double[] cthetaBuffer = ctheta.getDataXYAsDouble(0);
		IcyBufferedImage stheta = new IcyBufferedImage(_w, _h, 1, DataType.DOUBLE);
		double[] sthetaBuffer = stheta.getDataXYAsDouble(0);

		//7. Initialize pupil amplitude to one within bandwidth and phase to zero
		for(int x = 0; x < _w; x++)
		{
			for(int y = 0; y < _h; y++)
			{   
				double kxy = Math.sqrt( Math.pow((x-wc), 2) + Math.pow((y-hc), 2) );

				pupilMagBuffer[x + y * _w] = ((kxy < _kMax) ? 1 : 0); //Pupil bandwidth constraints
				pupilPhBuffer[x + y * _w] = 0; //Zero phase 
				sthetaBuffer[x + y * _w] = Math.sin( kxy * kSampling / _kObj );
				sthetaBuffer[x + y * _w] = (sthetaBuffer[x + y * _w]< 0) ? 0: sthetaBuffer[x + y * _w];
				cthetaBuffer[x + y * _w] = Double.MIN_VALUE + Math.sqrt(1 - Math.pow(sthetaBuffer[x + y * _w], 2));				
			}
		}
		//pupilImage.dataChanged();
		stheta.dataChanged();
		ctheta.dataChanged();		

		//8. Filter the pupil for antialiasing
		double[] gaussianKernel = Kernels1D.CUSTOM_GAUSSIAN.createGaussianKernel1D(_sigma).getData();
		double[][] tempPupil = new double[][]{ pupilMagBuffer };
		Convolution1D.convolve(tempPupil, _w, _h, gaussianKernel, gaussianKernel, null);	
		System.arraycopy(tempPupil[0], 0, pupilMagBuffer, 0, _w*_h);
		
		Sequence tpupil = new Sequence();
		tpupil.addImage(pupilImage);
		tpupil.setName("Estimated pupil");
		tpupil.setChannelName(0, "Magnitude");
		tpupil.setChannelName(1, "Phase");
		addSequence(tpupil);

		//KernelJAI gKernelJAI = new KernelJAI(_w, _h, Array2DUtil.doubleArrayToFloatArray(Kernels2D.CUSTOM.createCustomKernel2D(gaussianKernel, _w, _h, true)));
		//final int renderedop = ConvolveDescriptor.create(pupilImage, kernel, null);		

		//9. Iteration
		for(int n = 0; n<_nIter; n++)
		{
			IcyBufferedImage avgPupil = new IcyBufferedImage(_w, _h, 2, DataType.DOUBLE);
			double[] avgPupilReBuffer = avgPupil.getDataXYAsDouble(0); //Real
			double[] avgPupilImBuffer = avgPupil.getDataXYAsDouble(1); //imaginary

			for(int iz=0;iz<NSECTIONS;iz++)
			{
				//9a. Calculated Defocused pupil
				for(int x = 0; x < _w; x++)
				{
					for(int y = 0; y < _h; y++)
					{ 
						dpupilReBuffer[x + y * _w] = pupilMagBuffer[x + y * _w] * Math.cos(pupilPhBuffer[x + y * _w] + (defocus[iz] * _kObj * cthetaBuffer[x + y * _w]));
						dpupilImBuffer[x + y * _w] = pupilMagBuffer[x + y * _w] * Math.sin(pupilPhBuffer[x + y * _w] + (defocus[iz] * _kObj * cthetaBuffer[x + y * _w]));

					}
				}
				
				double[] psf2d = dpupilImage.getDataCopyCXYAsDouble();
				fftOp.complexForward(psf2d);


				//9b. Swap quadrants of PSF and update
				IcyBufferedImage psfCentered = new IcyBufferedImage(_w, _h, 2, DataType.DOUBLE);
				double[] psfReBuffer = psfCentered.getDataXYAsDouble(0);//Real
				double[] psfImBuffer = psfCentered.getDataXYAsDouble(1);//imaginary

				for(int x = 0; x < (wc+1); x++)
				{
					for(int y = 0; y < (hc+1); y++)
					{
						double r = Math.sqrt( Math.pow(x-wc, 2) + Math.pow(y-hc, 2) );
						psfReBuffer[x + y * _w] = psf2d[(((wc-x) + (hc-y) * _w)*2) + 0];	
						psfImBuffer[x + y * _w] = psf2d[(((wc-x) + (hc-y) * _w)*2) + 1];
						double psf = Double.MIN_VALUE + Math.pow(psfReBuffer[x + y * _w], 2) + Math.pow(psfImBuffer[x + y * _w], 2);
						//Update 
						psfReBuffer[x + y * _w] = ((r < _rMax) ? 1 : 0) * psfReBuffer[x + y * _w]  * (1 - _alpha - (_alpha * bgRemovedArray[cPlane + selectedPlanes[iz]][x + y * _w]/psf) ); 
						psfImBuffer[x + y * _w] = ((r < _rMax) ? 1 : 0) * psfImBuffer[x + y * _w]  * (1 + _alpha - (_alpha * bgRemovedArray[cPlane + selectedPlanes[iz]][x + y * _w]/psf) );

					}
					for(int y = hc+1; y < _h; y++)
					{
						double r = Math.sqrt( Math.pow(x-wc, 2) + Math.pow(y-hc, 2) );
						psfReBuffer[x + y * _w] = psf2d[(((wc-x) + (_h+ hc-y) * _w)*2) + 0];	
						psfImBuffer[x + y * _w] = psf2d[(((wc-x) + (_h+ hc-y) * _w)*2) + 1];
						double psf = Double.MIN_VALUE + Math.pow(psfReBuffer[x + y * _w], 2) + Math.pow(psfImBuffer[x + y * _w], 2);
						// Update 
						psfReBuffer[x + y * _w] = ((r < _rMax) ? 1 : 0) * psfReBuffer[x + y * _w]  * (1 - _alpha - (_alpha * bgRemovedArray[cPlane + selectedPlanes[iz]][x + y * _w]/psf) ); 
						psfImBuffer[x + y * _w] = ((r < _rMax) ? 1 : 0) * psfImBuffer[x + y * _w]  * (1 + _alpha - (_alpha * bgRemovedArray[cPlane + selectedPlanes[iz]][x + y * _w]/psf) );

					}

				}
				for(int x = (wc+1); x < _w; x++)
				{
					for(int y = 0; y < (hc+1); y++)
					{
						double r = Math.sqrt( Math.pow(x-wc, 2) + Math.pow(y-hc, 2) );
						psfReBuffer[x + y * _w] = psf2d[(((_w+wc-x) + (hc-y) * _w)*2) + 0];	
						psfImBuffer[x + y * _w] = psf2d[(((_w+wc-x) + (hc-y) * _w)*2) + 1];
						double psf = Double.MIN_VALUE + Math.pow(psfReBuffer[x + y * _w], 2) + Math.pow(psfImBuffer[x + y * _w], 2);
						// Update 
						psfReBuffer[x + y * _w] = ((r < _rMax) ? 1 : 0) * psfReBuffer[x + y * _w]  * (1 - _alpha - (_alpha * bgRemovedArray[cPlane + selectedPlanes[iz]][x + y * _w]/psf) ); 
						psfImBuffer[x + y * _w] = ((r < _rMax) ? 1 : 0) * psfImBuffer[x + y * _w]  * (1 + _alpha - (_alpha * bgRemovedArray[cPlane + selectedPlanes[iz]][x + y * _w]/psf) );
					}
					for(int y = hc+1; y < _h; y++)
					{
						double r = Math.sqrt( Math.pow(x-wc, 2) + Math.pow(y-hc, 2) );
						psfReBuffer[x + y * _w] = psf2d[(((_w+wc-x) + (_h+ hc-y) * _w)*2) + 0];	
						psfImBuffer[x + y * _w] = psf2d[(((_w+wc-x) + (_h+ hc-y) * _w)*2) + 1];
						double psf = Double.MIN_VALUE + Math.pow(psfReBuffer[x + y * _w], 2) + Math.pow(psfImBuffer[x + y * _w], 2);
						// Update 
						psfReBuffer[x + y * _w] = ((r < _rMax) ? 1 : 0) * psfReBuffer[x + y * _w]  * (1 - _alpha - (_alpha * bgRemovedArray[cPlane + selectedPlanes[iz]][x + y * _w]/psf) ); 
						psfImBuffer[x + y * _w] = ((r < _rMax) ? 1 : 0) * psfImBuffer[x + y * _w]  * (1 + _alpha - (_alpha * bgRemovedArray[cPlane + selectedPlanes[iz]][x + y * _w]/psf) );

					}
				}
				psfCentered.dataChanged();

				//9c. Return centered result
				psf2d = psfCentered.getDataCopyCXYAsDouble();

				// 9c. Calculate the experimentally updated pupil function
				double[] pupilArray = psf2d;
				fftOp.complexInverse(pupilArray, false);

				//9d. Correct for defocus and average the pupils
				for(int x = 0; x < _w; x++)
				{
					for(int y = 0; y < _h; y++)
					{ 
						dpupilReBuffer[x + y * _w] = Math.sqrt(Math.pow(pupilArray[((x + y * _w) * 2) + 0], 2) + Math.pow(pupilArray[((x + y * _w) * 2) + 1], 2)) * Math.cos(Math.atan2(pupilArray[((x + y * _w) * 2) + 1], pupilArray[((x + y * _w) * 2) + 0])-(defocus[iz] * _kObj * cthetaBuffer[x + y * _w]));
						dpupilImBuffer[x + y * _w] = Math.sqrt(Math.pow(pupilArray[((x + y * _w) * 2) + 0], 2) + Math.pow(pupilArray[((x + y * _w) * 2) + 1], 2)) * Math.sin(Math.atan2(pupilArray[((x + y * _w) * 2) + 1], pupilArray[((x + y * _w) * 2) + 0])-(defocus[iz] * _kObj * cthetaBuffer[x + y * _w]));

						avgPupilReBuffer[x + y * _w] = avgPupilReBuffer[x + y * _w] + dpupilReBuffer[x + y * _w];
						avgPupilImBuffer[x + y * _w] = avgPupilImBuffer[x + y * _w] + dpupilImBuffer[x + y * _w];

					}
				}
				dpupilImage.dataChanged();

			}
			ArrayMath.divide(avgPupilReBuffer, NSECTIONS);
			ArrayMath.divide(avgPupilImBuffer, NSECTIONS);			
			avgPupil.dataChanged();
			
			for(int x = 0; x < _w; x++)
			{
				for(int y = 0; y < _h; y++)
				{ 
					pupilMagBuffer[x + y * _w] =  Math.sqrt(Math.pow(avgPupilReBuffer[x + y * _w], 2) + Math.pow(avgPupilImBuffer[x + y * _w], 2));
					pupilPhBuffer[x + y * _w] = Math.atan2(avgPupilImBuffer[x + y * _w], avgPupilReBuffer[x + y * _w]);
				}
			}
			if(Math.IEEEremainder(n, 2) == 0)
			{
				tempPupil = new double[][]{ pupilMagBuffer };
				Convolution1D.convolve(tempPupil, _w, _h, gaussianKernel, gaussianKernel, null);	
				System.arraycopy(tempPupil[0], 0, pupilMagBuffer, 0, _w*_h);
				tempPupil = new double[][]{ pupilPhBuffer };
				Convolution1D.convolve(tempPupil, _w, _h, gaussianKernel, gaussianKernel, null);	
				System.arraycopy(tempPupil[0], 0, pupilPhBuffer, 0, _w*_h);
			}
			for(int x = 0; x < _w; x++)
			{
				for(int y = 0; y < _h; y++)
				{ 
					double kxy = Math.sqrt( Math.pow((x-wc), 2) + Math.pow((y-hc), 2) );
					pupilMagBuffer[x + y * _w] = ((kxy < _kMax) ? 1 : 0) * pupilMagBuffer[x + y * _w];
					pupilPhBuffer[x + y * _w] = ((kxy < _kMax) ? 1 : 0) * pupilPhBuffer[x + y * _w];
				}
			}
			MathUtil.normalize(pupilMagBuffer);
			pupilImage.dataChanged();
		}
		pupil.addImage(pupilImage);
		return pupil;
	}










	@Override
	public void clean() {
		// TODO Auto-generated by Icy4Eclipse
	}
}
