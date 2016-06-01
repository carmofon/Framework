/*
 * File created originally on Jun 14, 2005
 */
package nl.tudelft.bt.model.timeconstraint;

import java.io.Serializable;

import nl.tudelft.bt.model.Model;

/**
 * Stores information relative to timestep constraints
 * 
 * @author jxavier
 */
public class TimeStepConstraint implements Serializable {
	private long _realtime;
	private float _timeStep;
	private float _time;
	private int _iteration;
	private String _name;

	/**
	 * Initializes a new instance of TimeStepConstraint
	 * 
	 * @param t
	 *            the time step
	 * @param n
	 *            the name of the contraint
	 */
	public TimeStepConstraint(float t, String n) {
		_realtime = System.currentTimeMillis();
		_timeStep = t;
		_time = Model.model().getTime();
		_iteration = Model.model().getIterationCounter();
		_name = n;
	}

	/**
	 * Initializes a new instance of TimeStepConstraint with default value of
	 * _timeStep = infinity
	 */
	public TimeStepConstraint() {
		this(Float.POSITIVE_INFINITY, "Default");
	}

	/**
	 * Get a header for a output file
	 * 
	 * @return header for output file
	 */
	public static String getHeaderWithString() {
		return firstPartOfHeader() + "\t" + "name" + "\n";
	}

	/**
	 * @return
	 */
	public static String getHeaderWithoutString() {
		return firstPartOfHeader() + "\n";
	}

	/**
	 * @return
	 */
	private static String firstPartOfHeader() {
		return "iteration" + "\t" + "time" + "\t" + "timeStep" + "\t"
				+ "realTime";
	}
	/**
	 * @return
	 */
	public String writeWithName() {
		return writeFirstPart() + "\t" + _name + "\n";
	}
	/**
	 * @return
	 */
	public String writeWithoutName() {
		return writeFirstPart() + "\n";
	}
	/**
	 * @return
	 */
	private String writeFirstPart() {
		return _iteration + "\t" + _time + "\t" + _timeStep + "\t" + _realtime;
	}
	/**
	 * @param t1
	 * @param t2
	 * @return the instance with minimum time step
	 */
	public static TimeStepConstraint getMinimum(TimeStepConstraint t1,
			TimeStepConstraint t2) {
		if (t1._timeStep < t2._timeStep) {
			return t1;
		}
		return t2;
	}

	/**
	 * @param t1
	 * @param t2
	 * @param t3
	 * @return the instance with minimum time step
	 */
	public static TimeStepConstraint getMinimum(TimeStepConstraint t1,
			TimeStepConstraint t2, TimeStepConstraint t3) {
		return getMinimum(getMinimum(t1, t2), t3);
	}

	/**
	 * @return the value of the time step
	 */
	public float getTimeStep() {
		return _timeStep;
	}
	/**
	 * @param _name
	 *            The _name to set.
	 */
	public void setName(String _name) {
		this._name = _name;
	}
	/**
	 * @param step
	 *            The _timeStep to set.
	 */
	public void setTimeStep(float step) {
		_timeStep = step;
	}
}