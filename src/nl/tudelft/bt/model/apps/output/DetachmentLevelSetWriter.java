package nl.tudelft.bt.model.apps.output;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.*;
/**
 * Writes the detachment level set matrix before actual detachment is carried
 * out at each iteration. The level set matrix is stored as a matrix which
 * represents the detachment times of each entry
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class DetachmentLevelSetWriter extends StateWriter {
	static final private String DETACHLSDIR = "detachmentLevelSet";
	/*
	 * (non-Javadoc)
	 */
	public void write() throws ModelException {
		try {
			Model.model().writeDetachmentLevelSet(
					confirmSubDirectoryExists(DETACHLSDIR));
		} catch (Exception e) {
			throw new ModelIOException("Error trying to write"
					+ " detachment level set");
		}
	}
}