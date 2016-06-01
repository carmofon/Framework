package nl.tudelft.bt.model.work.tumordemonstration;
import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.output.StateWriter;
import nl.tudelft.bt.model.exceptions.ModelException;
import nl.tudelft.bt.model.exceptions.ModelIOException;
import nl.tudelft.bt.model.povray.Povray3DScene;
import nl.tudelft.bt.model.util.ZipArchive;
/**
 * Writer for povray files (for 3D rendering of scenes using povray software)
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class PovRayShovingWriter extends StateWriter {
	private static PovRayShovingWriter _instance;
	
	private int iterationNumber = 0;
	private int previousIterationNumber = -1;
	private int shovingFrame = 0;
	
	static final private String POVRAYDIR = "povray";
	private Povray3DScene _povRay;
	/*
	 * (non-Javadoc)
	 */
	public void write() throws ModelException {
		try {
			if (Model.model().getIterationCounter() > previousIterationNumber) {
				previousIterationNumber = Model.model().getIterationCounter();
				shovingFrame = 0;
			} else
				shovingFrame++;
			
			//contruct directory name
			String dir = confirmSubDirectoryExists(POVRAYDIR);
			// if instance of PovRay is not initialize, initialize now and
			// create header and footer files and add to zip archive
			if (_povRay == null) {
				_povRay = new Povray3DScene();
			}
			// check if directory is empty, and if so create the header
			// and footer
			File directory = new File(dir);
			String[] ls = directory.list();
			if (ls.length < 1)
				_povRay.writePovrayIncFiles(dir);
			// write the state to file (just particles, using the include files
			String fn = dir + "it" + getFormattedIterationNumber()
					+ "_" + Model.model().getFormatedNumber(shovingFrame) + ".pov";
			File f = _povRay.writeModelStateToPovrayFile(fn);
			ZipArchive zipArchive = new ZipArchive(fn + ".zip"); 
			zipArchive.addToZipArchiveAndDelete(f);
			zipArchive.close();
		} catch (IOException e) {
			throw new ModelIOException("Error trying to write povRayFile");
		}
	}
	
	public static void setInstance(PovRayShovingWriter s) {
		_instance = s;
	}
	
	public static void writeNow() {
		try {
			_instance.write();
		} catch (ModelException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}	
}
