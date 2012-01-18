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

public class MacroscopePSF extends EzPlug {
	private EzVarInteger _w;
	private EzVarInteger _h;
	private EzVarInteger _z;
	private EzVarDouble _xySampling;
	private EzVarDouble _zSampling;
	private EzVarDouble _indexImmersion;
	private EzVarDouble _objNA;
	private EzVarDouble _zoomNA;
	private EzVarInteger _lem;	
	private EzVarInteger _xofactor;
	private EzVarInteger _yofactor;
	

	@Override
	protected void initialize() {
		_w = new EzVarInteger("Image width in pixels");				
		_h = new EzVarInteger("Image height in pixels");
		_z = new EzVarInteger("Number of slices in the volume");
		_xySampling = new EzVarDouble("Image pixel spacing, in nm");
		_zSampling = new EzVarDouble("Slice spacing (z), in nm");
		_indexImmersion = new EzVarDouble("Refractive index of the medium between lens and cover slip (default air)");
		_objNA = new EzVarDouble("Effective numerical aperture of the objective lens");
		_zoomNA = new EzVarDouble("Effective numerical aperture of the zoom system");
		_lem = new EzVarInteger("Emission peak wavelength, in nm");
		_xofactor = new EzVarInteger("Relative x displacement (in pixels)", 0, 200, 1);
		_yofactor = new EzVarInteger("Relative y displacement (in pixels)", 0, 200, 1);
		        
		// Set the default values
        _w.setValue(MacroCalculator.DEFAULT_W);
        _h.setValue(MacroCalculator.DEFAULT_H);
        _z.setValue(MacroCalculator.DEFAULT_Z);
        _xySampling.setValue(MacroCalculator.DEFAULT_XYSAMPLING);
        _zSampling.setValue(MacroCalculator.DEFAULT_ZSAMPLING);
        _indexImmersion.setValue(MacroCalculator.DEFAULT_INDEXIMMERSION);
        _objNA.setValue(MacroCalculator.DEFAULT_OBJNA);
        _zoomNA.setValue(MacroCalculator.DEFAULT_ZOOMNA);
        _lem.setValue(MacroCalculator.DEFAULT_LEM);
        _xofactor.setValue(MacroCalculator.DEFAULT_XOFACTOR);
        _yofactor.setValue(MacroCalculator.DEFAULT_YOFACTOR);
        EzGroup parameterGroup = new EzGroup("Enter MACROscope settings", _w, _h, _z, _xySampling, _zSampling, _indexImmersion, _objNA, _zoomNA, _lem, _xofactor, _yofactor);
		addEzComponent(parameterGroup);        
	}

	@Override
	protected void execute() {
		MacroCalculator parameters = new MacroCalculator();
		parameters.setW(_h.getValue());
		parameters.setH(_h.getValue());
		parameters.setZ(_z.getValue());
		parameters.setXYSAMPLING(_xySampling.getValue());
		parameters.setZSAMPLING(_zSampling.getValue());
		parameters.setIndexImmersion(_indexImmersion.getValue());
		parameters.setOBJNA(_objNA.getValue());
		parameters.setZOOMNA(_zoomNA.getValue());
		parameters.setLEM(_lem.getValue());
		parameters.setXOFACTOR(_xofactor.getValue());
		parameters.setYOFACTOR(_yofactor.getValue());
		addSequence(parameters.compute());		// TODO Auto-generated by Icy4Eclipse
		//MessageDialog.showDialog("test is working fine !");
	}

	@Override
	public void clean() {
		// TODO Auto-generated by Icy4Eclipse
	}
}
