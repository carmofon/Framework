package nl.tudelft.bt.model.apps.output;

import java.io.File;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.profiles1d.Profile;

/**
 * Writes a profile to disk
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ProfileWriter extends StateWriter {
	static final private String PROFILEDIR = "profiles";
	private Profile _profile;

	/**
	 * @param profile
	 * @throws SystemEditViolationException
	 */
	public ProfileWriter(Profile profile) throws SystemEditViolationException {
		super();
		_profile = profile;
	}

	/*
	 * (non-Javadoc)
	 */
	public void write() throws ModelException {
		String baseDir = confirmSubDirectoryExists(PROFILEDIR);
		String subDir = "iteration"
				+ Model.model().getFormatedIterationNumber();
		// create directory for this iteration
		File dir = new java.io.File(baseDir + "/" + subDir);
		try {
			dir.mkdir();
		} catch (SecurityException e) {
			throw new SystemEditViolationException("directory " + dir
					+ " could not be created for security reasons");
		}
		//create file
		String fileName = PROFILEDIR + "/" + subDir + "/" + _profile.getName()
				+ ".txt";
		//write to file
		appendToFile(fileName, _profile.getFormatedTable());

	}
}