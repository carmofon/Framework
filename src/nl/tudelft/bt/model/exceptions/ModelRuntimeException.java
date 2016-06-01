/*
 * Created on 6-feb-2004 by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.exceptions;
/**
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ModelRuntimeException extends RuntimeException {
	/**
	 *  
	 */
	public ModelRuntimeException() {
		super();
	}
	/**
	 * @param message
	 */
	public ModelRuntimeException(String message) {
		super(message);
	}
	/**
	 * Create an exception to be thrown when a given number is not valid in
	 * runtime
	 * 
	 * @param f
	 */
	public ModelRuntimeException(float f) {
		super(f + " is not a valid number");
	}
}
