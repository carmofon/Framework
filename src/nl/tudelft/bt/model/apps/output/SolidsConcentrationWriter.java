package nl.tudelft.bt.model.apps.output;


import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.*;

/**
 * Writes concentration fields of all solids (bacteria species and other
 * shovable components) to file in disk
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class SolidsConcentrationWriter extends StateWriter {
	static final public String SOLIDDIR = "solids";

	/*
	 * (non-Javadoc)
	 */
	public void write() throws ModelException {
		try {
			Model.model().writeSolidConcentrationsToFile(
					confirmSubDirectoryExists(SOLIDDIR));
		} catch (Exception e) {
			throw new ModelIOException("Error trying to write"
					+ " solid concentrations");
		}

	}
}