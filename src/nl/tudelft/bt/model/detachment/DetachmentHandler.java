package nl.tudelft.bt.model.detachment;

public interface DetachmentHandler {

	/**
	 * Method that determines if detachment is off, allowing the simulation to
	 * avoid computing the detachment routines, which results in improved
	 * performance.
	 * 
	 * 
	 * @return true if detachment is not working
	 */
	public abstract boolean detachmentIsOff();

}