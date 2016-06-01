package nl.tudelft.bt.model.apps.output;

import java.io.File;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.exceptions.ModelException;
import nl.tudelft.bt.model.exceptions.ModelIOException;
import nl.tudelft.bt.model.exceptions.SystemEditViolationException;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.povray.Povray3DScene;
import nl.tudelft.bt.model.util.ZipArchive;

/**
 * Writes image directly to file
 * 
 * @author Joao Xavier (xavierj@mskcc.org)
 */
public class ImageWriter extends StateWriter {
	static final private String FRAMESDIR = "frames";
	private ImageMaker _imMaker;
	private String _dir;
	private String _stringToAppend ;

	public ImageWriter(SoluteSpecies s) throws SystemEditViolationException {
		super();
		_dir = confirmSubDirectoryExists(FRAMESDIR);
		try {
			_stringToAppend = s.getName();
			_imMaker = new ImageMaker(s, _dir + _stringToAppend + ".mov");
		} catch (IOException e) {
			throw new ModelIOException("Error trying to write image");
		}
	}

	public ImageWriter() throws SystemEditViolationException {
		super();
		_dir = confirmSubDirectoryExists(FRAMESDIR);
		try {
			_stringToAppend = "sim";
			_imMaker = new ImageMaker(null, _dir + _stringToAppend + ".mov");
		} catch (IOException e) {
			throw new ModelIOException("Error trying to write image");
		}
	}

	/*
	 * (non-Javadoc)
	 */
	public void write() throws ModelException {
		try {
			// write the state to file (just particles, using the include files
			String fn = _dir + _stringToAppend + getFormattedIterationNumber();
			_imMaker.writeFrame(fn);
		} catch (IOException e) {
			throw new ModelIOException("Error trying to write image");
		}
	}
}
