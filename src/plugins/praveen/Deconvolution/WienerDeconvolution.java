package plugins.praveen.Deconvolution;

import icy.sequence.Sequence;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.ezplug.EzVarText;

public class WienerDeconvolution extends EzPlug{
	/**
	 *	Wiener Deconvolution
	 *
	 * @author Praveen Pankajakshan
	 * 
	 * @see "Digital Image Processing", R. C. Gonzalez & R. E. Woods, Addison-Wesley Publishing Company, Inc., 1992.
	 * 
	 */
	EzVarSequence observed = new EzVarSequence("Observed stack");
	EzVarSequence psf = new EzVarSequence("PSF");
	EzVarDouble nsr = new EzVarDouble("Noise-to-signal ratio", 0, 0, 100, 0.001);
	EzVarText dim = new EzVarText("Dimension", new String[] { "2-D", "3-D" },0, false);
	EzVarDouble ratio = new EzVarDouble("Z/X rezolution Ratio", 1, 0.001, 100, 0.001);
	
	@Override
	protected void initialize() {
		// TODO Auto-generated method stub
		addEzComponent(observed);
		addEzComponent(psf);
		addEzComponent(dim);
		addEzComponent(nsr);
		addEzComponent(ratio);

	}
	
	@Override
	protected void execute() {
		Sequence deconvolved = null;
		if(dim.getValue().equalsIgnoreCase("2-D"))
		{
			deconvolved = Wiener2D(observed.getValue(), psf.getValue(), nsr.getValue(), ratio.getValue());
		}
		else
		{
			deconvolved = Wiener3D(observed.getValue(), psf.getValue(), nsr.getValue(), ratio.getValue());
		}
		addSequence(deconvolved);


	}
	
	private Sequence Wiener3D(Sequence inSeq, Sequence h, double ratioNoiseSignal,
			double ratio) {
		int nxSeq = inSeq.getSizeX();
		int nySeq = inSeq.getSizeY();
		int nzSeq = inSeq.getSizeZ();
		
		int nxPSF = h.getSizeX();
		int nyPSF = h.getSizeY();
		int nzPSF = h.getSizeZ();
		
		return outSeq;
	}

	private Sequence Wiener2D(Sequence value, Sequence value2, Double value3,
			Double value4) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void clean() {
		// TODO Auto-generated method stub

	}

	

	

}
