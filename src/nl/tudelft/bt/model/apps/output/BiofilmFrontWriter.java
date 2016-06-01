package nl.tudelft.bt.model.apps.output;
import java.io.IOException;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.*;
/**
 * Writes a boolean matrix with the biofilm front.
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class BiofilmFrontWriter extends StateWriter {
	static final private String DETACHLSDIR = "biofilmFront";
	/*
	 * (non-Javadoc)
	 */
	public void write() throws ModelException {
		try {
			Model.model().writeBiofilmFront(
					confirmSubDirectoryExists(DETACHLSDIR));
		} catch (IOException e) {
			throw new ModelIOException("Error trying to write"
					+ " biofilm front");
		}
	}
}