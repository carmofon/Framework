package nl.tudelft.bt.model.exceptions;

/**
 * General excpetion class for the multidimensional biofilm model
 * 
 * @author João Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ModelException extends Exception {

	/**
	 * 
	 */
	public ModelException() {
		super();
	}

	/**
	 * @param arg0
	 */
	public ModelException(String arg0) {
		super(arg0);
	}
}
