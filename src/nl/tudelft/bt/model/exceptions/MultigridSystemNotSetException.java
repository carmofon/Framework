package nl.tudelft.bt.model.exceptions;

/**
 * Exception to be thrown when an application is not in editable
 * mode and changes to parameters are requested 
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class MultigridSystemNotSetException extends ModelException {
	/**
	 * 
	 */
	public MultigridSystemNotSetException() {
		super("Tried creation of MG variable but MG sytem not defined yet");
	}
}
