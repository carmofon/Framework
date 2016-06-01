package nl.tudelft.bt.model.exceptions;

/**
 * Exception to be thrown in case of inconsistent number
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class NonMatchingNumberException extends ModelRuntimeException {

	/**
	 * 
	 */
	public NonMatchingNumberException() {
		super();
	}

	/**
	 * @param arg0
	 */
	public NonMatchingNumberException(String arg0) {
		super(arg0);
	}
}
