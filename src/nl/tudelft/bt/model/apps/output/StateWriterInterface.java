package nl.tudelft.bt.model.apps.output;

import nl.tudelft.bt.model.exceptions.ModelException;

/**
 * Implements an interface for the sate wirters, allowing the decorator
 * patter to be applied (i.e. allowing creation of TimedStateWrites, etc.)
 */
public interface StateWriterInterface {
	/**
	 * Write a set of model ooutput to a given directory
	 */
	public abstract void write() throws ModelException;
	/**
	 * Force writting to disk, for debugging purposes
	 */
	public abstract void dump();
}