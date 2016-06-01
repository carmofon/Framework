package nl.tudelft.bt.model.exceptions;

/**
 * Exception to be thrown when an application is not in editable
 * mode and changes to parameters are requested 
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class SystemEditViolationException extends ModelRuntimeException{
	/**
	 * 
	 */
	public SystemEditViolationException(String msg) {
		super(msg);
	}
}
