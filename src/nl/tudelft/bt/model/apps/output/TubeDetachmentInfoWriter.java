package nl.tudelft.bt.model.apps.output;


import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.*;

/**
 * Writes concentration fields of all solutes to file in disk
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class TubeDetachmentInfoWriter extends StateWriter {
	static final private String SUBDIR = "tube-detachment";

	/*
	 * (non-Javadoc)
	 */
	public void write() throws ModelException {
		try {
			Model.model().writeTubeDetachmentAlgorithmOutput(
					confirmSubDirectoryExists(SUBDIR));
		} catch (Exception e) {
			e.printStackTrace();
			throw new ModelIOException("Error trying to write to "
					+ SUBDIR);
		}

	}
}