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
import plugins.adufour.ezplug.EzVarText;
//import plugins.adufour.ezplug.EzVarSequence;

public class ConfocalPSF extends EzPlug {
	private EzVarInteger _w;
	private EzVarInteger _h;
	private EzVarInteger _z;
	private EzVarText _mName;
	private EzVarDouble _indexImmersion;
	private EzVarDouble _na;
	private EzVarInteger _mObj;
	private EzVarInteger _lex;
	private EzVarInteger _lem;	
	private EzVarDouble _indexSpecimen;
	private EzVarDouble _xySampling;
	private EzVarDouble _zSampling;
	private EzVarDouble _pSize;	
	private EzVarDouble _depth;	

	@Override
	protected void initialize() {		
		_mName = new EzVarText("Choose your microscope", new String[] { "Biorad MRC 500/600/1024", "Biorad Radiance", "Leica TCS 4D/SP1/NT", "Leica SP2", "Leica SP5", "TE2000-E C1 Head", "Ti-E Perfect Focus A1R", "Olympus FV10i", "Olympus FV300/FVX", "Olympus FV500", "Olympus FV1000", "Visitech Infinity", "Zeiss LSM410", "Zeiss LSM510", "Zeiss LSM700", "Zeiss LSM710" }, 0, false);
		_na = new EzVarDouble("Actual numerical aperture of objective");
		_mObj = new EzVarInteger("Objective Magnification");
		_indexImmersion = new EzVarDouble("Refractive index of the lens immersion medium (default oil)");
		_pSize = new EzVarDouble("Pinhole radius in nm");//http://www.svi.nl/BackprojectedPinholeCalculator
		_lex = new EzVarInteger("Excitation peak wavelength, in nm");
		_lem = new EzVarInteger("Emission peak wavelength, in nm");
		_indexSpecimen = new EzVarDouble("Refractive index of the specimen mounting medium (default water)");	
		_xySampling = new EzVarDouble("Image pixel spacing, in nm");
		_zSampling = new EzVarDouble("Slice spacing (z), in nm");	
		_depth = new EzVarDouble("Depth of imaging (under coverslip) in nm");
		_w = new EzVarInteger("Image width in pixels");				
		_h = new EzVarInteger("Image height in pixels");
		_z = new EzVarInteger("Number of slices in volume");
		
		// Set the default values
		_mName.setValue(ConfocalCalculator.DEFAULT_MNAME);
		_na.setValue(ConfocalCalculator.DEFAULT_NA);
		_mObj.setValue(ConfocalCalculator.DEFAULT_MOBJ);
        _indexImmersion.setValue(ConfocalCalculator.DEFAULT_INDEXIMMERSION);
        _pSize.setValue(ConfocalCalculator.DEFAULT_PSIZE);
        _lex.setValue(ConfocalCalculator.DEFAULT_LEX);
        _lem.setValue(ConfocalCalculator.DEFAULT_LEM);
        _indexSpecimen.setValue(ConfocalCalculator.DEFAULT_INDEXSP);
        _xySampling.setValue(ConfocalCalculator.DEFAULT_XYSAMPLING);
        _zSampling.setValue(ConfocalCalculator.DEFAULT_ZSAMPLING);
        _depth.setValue(ConfocalCalculator.DEFAULT_DEPTH);
        _w.setValue(ConfocalCalculator.DEFAULT_W);
        _h.setValue(ConfocalCalculator.DEFAULT_H);
        _z.setValue(ConfocalCalculator.DEFAULT_Z);        
        EzGroup parameterGroup = new EzGroup("Enter Confocal Microscope settings", _mName, _na, _mObj, _indexImmersion, _pSize, _lex, _lem, _indexSpecimen, _xySampling, _zSampling, _depth, _w, _h, _z);
		addEzComponent(parameterGroup);        
	}

	@Override
	protected void execute() {
		ConfocalCalculator parameters = new ConfocalCalculator();
		parameters.setMNAME(_mName.getValue());
		parameters.setNA(_na.getValue());
		parameters.setMOBJ(_mObj.getValue());
		parameters.setIndexImmersion(_indexImmersion.getValue());
		parameters.setPSIZE(_pSize.getValue());
		parameters.setLEM(_lem.getValue());
		parameters.setLEX(_lex.getValue());
		parameters.setIndexSp(_indexSpecimen.getValue());
		parameters.setXYSAMPLING(_xySampling.getValue());
		parameters.setZSAMPLING(_zSampling.getValue());
		parameters.setDEPTH(_depth.getValue());
		parameters.setW(_h.getValue());
		parameters.setH(_h.getValue());
		parameters.setZ(_z.getValue());		
		addSequence(parameters.compute());		// TODO Auto-generated by Icy4Eclipse
		//MessageDialog.showDialog("test is working fine !");
	}

	@Override
	public void clean() {
		// TODO Auto-generated by Icy4Eclipse
	}
}