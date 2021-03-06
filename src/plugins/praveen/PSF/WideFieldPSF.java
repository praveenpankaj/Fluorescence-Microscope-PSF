package plugins.praveen.PSF;

import plugins.adufour.ezplug.EzPlug;
//import icy.gui.dialog.MessageDialog;
//import icy.image.IcyBufferedImage;
//import icy.type.TypeUtil;

//import icy.gui.frame.GenericFrame;
//import icy.gui.frame.progress.AnnounceFrame;
//import icy.sequence.Sequence;

//import java.awt.event.ActionEvent;
//import java.awt.event.ActionListener;

//import javax.swing.JTextPane;

//import plugins.adufour.ezplug.EzButton;
import plugins.adufour.ezplug.EzGroup;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarInteger;
//import plugins.adufour.ezplug.EzVarSequence;

public class WideFieldPSF extends EzPlug {
	private EzVarInteger _w;
	private EzVarInteger _h;
	private EzVarInteger _z;
	private EzVarDouble _indexImmersion;
	private EzVarDouble _na;
	private EzVarInteger _lem;	
	private EzVarDouble _indexSpecimen;
	private EzVarDouble _xySampling;
	private EzVarDouble _zSampling;
	private EzVarDouble _depth;

	@Override
	protected void initialize() {
		_w = new EzVarInteger("Image width in pixels");				
		_h = new EzVarInteger("Image height in pixels");
		_z = new EzVarInteger("Number of slices in volume");
		_xySampling = new EzVarDouble("Image pixel spacing, in nm");
		_zSampling = new EzVarDouble("Slice spacing (z), in nm");
		_indexImmersion = new EzVarDouble("Refractive index of the lens immersion medium (default oil)");
		_na = new EzVarDouble("Effective numerical aperture");
		_lem = new EzVarInteger("Emission peak wavelength, in nm");
		_indexSpecimen = new EzVarDouble("Refractive index of the specimen mounting medium (default water)");		
		_depth = new EzVarDouble("Depth of imaging (under coverslip) in nm");
		        
		// Set the default values
        _w.setValue(PSFCalculator.DEFAULT_W);
        _h.setValue(PSFCalculator.DEFAULT_H);
        _z.setValue(PSFCalculator.DEFAULT_Z);
        _indexImmersion.setValue(PSFCalculator.DEFAULT_INDEXIMMERSION);
        _na.setValue(PSFCalculator.DEFAULT_NA);
        _lem.setValue(PSFCalculator.DEFAULT_LEM);
        _indexSpecimen.setValue(PSFCalculator.DEFAULT_INDEXSP);
        _xySampling.setValue(PSFCalculator.DEFAULT_XYSAMPLING);
        _zSampling.setValue(PSFCalculator.DEFAULT_ZSAMPLING);
        _depth.setValue(PSFCalculator.DEFAULT_DEPTH);
        
        EzGroup parameterGroup = new EzGroup("Enter Widefield Microscope settings", _w, _h, _z, _xySampling, _zSampling, _indexImmersion, _na, _lem, _indexSpecimen, _depth);
		addEzComponent(parameterGroup);        
	}

	@Override
	protected void execute() {
		PSFCalculator parameters = new PSFCalculator();
		parameters.setW(_h.getValue());
		parameters.setH(_h.getValue());
		parameters.setZ(_z.getValue());
		parameters.setIndexImmersion(_indexImmersion.getValue());
		parameters.setNA(_na.getValue());
		parameters.setLEM(_lem.getValue());
		parameters.setIndexSp(_indexSpecimen.getValue());
		parameters.setXYSAMPLING(_xySampling.getValue());
		parameters.setZSAMPLING(_zSampling.getValue());
		parameters.setDEPTH(_depth.getValue());
		
		addSequence(parameters.compute());		// TODO Auto-generated by Icy4Eclipse
		//MessageDialog.showDialog("test is working fine !");
	}

	@Override
	public void clean() {
		// TODO Auto-generated by Icy4Eclipse
	}
}
