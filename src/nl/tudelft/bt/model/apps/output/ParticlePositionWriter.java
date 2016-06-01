package nl.tudelft.bt.model.apps.output;


import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.*;

/**
 * Writes the position and other properties (according to biomass particle
 * method toString) to file
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ParticlePositionWriter extends StateWriter {
	static final private String PARTICLEDIR = "particles";

	/*
	 * (non-Javadoc)
	 */
	public void write() throws ModelException {
		try {
			Model.model().writeParticlePositionsToFile(
					confirmSubDirectoryExists(PARTICLEDIR));
		} catch (Exception e) {
			throw new ModelIOException("Error trying to write"
					+ " particle concentrations");
		}

	}
}