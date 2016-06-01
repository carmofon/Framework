package nl.tudelft.bt.model.apps.output;


import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.*;

/**
 * Writes concentration fields of all solutes to file in disk
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class SoluteConcentrationWriter extends StateWriter {
	static final private String SOLUTEDIR = "solutes";

	/*
	 * (non-Javadoc)
	 */
	public void write() throws ModelException {
		try {
			Model.model().writeSoluteConcentrationsToFile(
					confirmSubDirectoryExists(SOLUTEDIR));
		} catch (Exception e) {
			throw new ModelIOException("Error trying to write"
					+ " solute concentrations");
		}

	}
}