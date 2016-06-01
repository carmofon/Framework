package nl.tudelft.bt.model.apps.output;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.TubeVelocityField;

/**
 * Writes the velocity and shear stress fields to disk
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class FlowWriter extends StateWriter {
	static final private String FLOWDIR = "flow";

	/*
	 * (non-Javadoc)
	 */
	public void write() throws ModelException {
		try {
			Model.model().writeMultigridVariableToFile(
					confirmSubDirectoryExists(FLOWDIR),
					TubeVelocityField.getInstance());
			Model.model().writeMultigridVariableToFile(
					confirmSubDirectoryExists(FLOWDIR),
					TubeVelocityField.getInstance().getShearStress());
		} catch (Exception e) {
			throw new ModelIOException("Error trying to write"
					+ " velocity fields");
		}
	}
}