package nl.tudelft.bt.model.work.will.cell;

import nl.tudelft.bt.model.ContinuousCoordinate;

/**
 * Interface for implementing a motile cell.
 * 
 * @author wkc
 *
 */

public interface migration {
	
	public float calculateMovementTime();
	
	public ContinuousCoordinate calculateDirection();
	
	public void migrate();
}
