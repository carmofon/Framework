/* 
 * Created on 10-feb-2004 
 * by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.detachment.levelset.functions;

import java.io.Serializable;

import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.detachment.DetachmentHandler;

/**
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public abstract class DetachmentSpeedFunction implements Serializable, DetachmentHandler {
	/**
	 * @param c
	 * @return the value of the detachment speed [ um h^-1]
	 */
	public abstract float getValue(ContinuousCoordinate c);

	/* (non-Javadoc)
	 * @see nl.tudelft.bt.model.detachment.levelset.functions.DetachmentHandler#detachmentIsOff()
	 */
	public abstract boolean detachmentIsOff();
}