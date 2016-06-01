/*
 * File created originally on Jul 21, 2005
 */
package nl.tudelft.bt.model.profiles1d;


/**
 * @author jxavier
 */
public interface Profile {
	/**
	 * @return a string with a table with the values in the profile
	 */
	public abstract String getFormatedTable();
	/**
	 * @return the name of the profile
	 */
	public abstract String getName();	
}