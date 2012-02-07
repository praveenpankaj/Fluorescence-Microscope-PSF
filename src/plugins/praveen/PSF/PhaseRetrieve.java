package plugins.praveen.PSF;

import icy.gui.dialog.MessageDialog;
import icy.gui.frame.progress.AnnounceFrame;
import icy.image.IcyBufferedImage;
import icy.sequence.Sequence;
import icy.type.DataType;
import icy.type.collection.array.Array1DUtil;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarFloat;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.filtering.Convolution1D;
import plugins.adufour.filtering.Kernels1D;
import plugins.adufour.projection.Projection;
import icy.math.ArrayMath;
import icy.math.MathUtil;
import javax.media.jai.BorderExtender;
import javax.media.jai.JAI;
import javax.media.jai.RenderedOp;
import javax.media.jai.operator.BorderDescriptor;
import javax.swing.SwingConstants;


public class PhaseRetrieve extends EzPlug {
	EzVarSequence _input = new EzVarSequence("Choose the 3D fluorescence bead image");
	EzVarFloat _xySampling = new EzVarFloat("Image pixel spacing, in nm", (float)92.00, (float)10.00, (float)50000.00, (float)0.1);
	EzVarFloat _zSampling = new EzVarFloat("Slice spacing (z), in nm", (float)277.00, (float)10.00, (float)50000.00, (float)0.1);
	EzVarFloat _objNA = new EzVarFloat("Effective numerical aperture of the objective lens", (float)1.4, (float)0.1, (float)4.00, (float)0.01);
	EzVarFloat _indexImmersion = new EzVarFloat("Refractive index of the lens immersion medium", (float)1.518, (float)1.00, (float)4.00, (float)0.01);
	EzVarInteger _lem = new EzVarInteger("Emission peak wavelength, in nm", 520, 405, 750, 1);
	EzVarInteger _bgd = new EzVarInteger("Mean background fluorescence intensity", 0, 0, 100000, 1);
	EzVarFloat _sigma = new EzVarFloat("Gaussian filter parameter", (float)0.5, (float)0.1, (float)1.00, (float)0.01);
	EzVarFloat _alpha = new EzVarFloat("Step size for the iterative algorithm", (float)0.6, (float)0.5, (float)1.00, (float)0.01); 
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
		addSequence(pupil);		
		pupil.setName("Estimated Back Aperture Pupil");
		pupil.setChannelName(0, "Magnitude");
		pupil.setChannelName(1, "Phase");
		//MessageDialog.showDialog("Test is working fine!");
	}

	public Sequence estimatepupil(Sequence sequence, float _xySampling, float _zSampling, float _objNA, float _indexImmersion, int _lem, int _bgd, float _sigma, float _alpha, int _nIter) {
		// TODO Auto-generated method stub
		Sequence pupil = new Sequence();
		Sequence resizedSeq = new Sequence();
		Sequence selectedSeq = new Sequence();

		int _w = sequence.getSizeX();
		int _h = sequence.getSizeY();
		int _z = sequence.getSizeZ();
		double[][] seqArray = new double[_z][_w*_h];
		double[][] bgRemovedArray = new double[_z][_w*_h];

		//1. Calculate the parameter necessary for the algorithm
		double _lambdaObj = _lem/_indexImmersion; //Wavelength inside the lens
		double _kObj = 2*Math.PI/_lambdaObj; //Wavenumber inside the lens
		double _kMax = 2*Math.PI*_objNA/_lambdaObj; //Maximum permissible frequency
		double _xyMax = 3*0.061*_lambdaObj/_objNA;//Maximum spread

		// 1. Resize the input data to make it a square image
		int leftPad = 0;
		int rightPad = 0;
		int topPad = 0;
		int botPad = 0;
		if(_w>_h)
		{ 
			int dh = _w-_h;
			_h = _w;
			if(Math.IEEEremainder(dh, 2)==0)
			{
				topPad = dh/2;
				botPad = dh/2;

			}
			else
			{
				topPad = (int) Math.ceil(dh/2);
				botPad = (int) Math.floor(dh/2);
			}
		}
		else
		{
			int dw = _h -_w;
			_w = _h;
			if(Math.IEEEremainder(dw, 2)==0)
			{
				leftPad = dw/2;
				rightPad = dw/2;

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
			resizedSeq.addImage(resizedImage);
		}
		double kSampling = 2*Math.PI/(_w*_xySampling);

		//2. Find central plane
		//Sequence zMaxProj = new Sequence();
		int cPlane = 0;
		double[] zMaxIntensity = new double[_z];
		double maxIntensity=0;
		Sequence zMaxProj = Projection.zProjection(sequence, Projection.ProjectionType.MAX, true);


		// 3. Remove Mean Background FLuorescence Intensity
		for(int iz = 0;iz<_z;iz++)
		{		
			IcyBufferedImage zMaxProjImage = zMaxProj.getImage(0, iz, 0);
			zMaxIntensity[iz] = zMaxProjImage.getComponentUserMaxValue(0);
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
					bgRemovedArray[iz][ix + iy*_h] = seqArray[iz][ix + iy*_h]-_bgd;
					bgRemovedArray[iz][ix + iy*_h] = ((bgRemovedArray[iz][ix + iy*_h] < 0) ? 0 : bgRemovedArray[iz][ix + iy*_h]);					
				}
			}			
		}
		//4. Display the focal plane information
		new AnnounceFrame("Detected focal plane at " + cPlane+1 + "th slice.");
		int[] selectedPlanes = new int[]{cPlane+1-15, cPlane+1-2, cPlane+1+2, cPlane+1+15};
		float[] defocus = new float[4];
		ArrayMath.multiply(ArrayMath.subtract(selectedPlanes, cPlane), _zSampling, defocus);	
		//5. Initialize Pupil Function
		for(int ix = 0;ix<_w;ix++)
		{
			for(int iy = 0;iy<_h;iy++)
			{
				pupilMagArray[0][ix + iy*_h] = 0;
				pupilPhArray[0][ix + iy*_h] = 0;
			}
		}



		double[] gaussianKernel = Kernels1D.CUSTOM_GAUSSIAN.createGaussianKernel1D(_sigma).getData();

		return null;
	}

	@Override
	public void clean() {
		// TODO Auto-generated by Icy4Eclipse
	}
}
