/* 
 * Created on 16-feb-2004 
 * by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.exceptions;

/**
 * Exception thrown when cells are not growing for lack of nutrients
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class GrowthProblemException extends ModelRuntimeException {

	/**
	 * 
	 */
	public GrowthProblemException() {
		super("Cells is not growing for lack of nutrients");
	}
	/**
	 * @param message
	 */
	public GrowthProblemException(String message) {
		super(message);
	}
}
