/*
 * Created on 6-feb-2004 by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.apps.output;
import java.io.Serializable;

import nl.tudelft.bt.model.util.UnlimitedFloatArray;
/**
 * Interface for varibles for plots
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public abstract class VariableSeries implements Serializable {
	private String _name;
	private String _xLabel;
	private String _yLabel;
	private UnlimitedFloatArray _x;
	private UnlimitedFloatArray _y;
	/**
	 *  
	 */
	public VariableSeries(String name, String xlabel, String ylabel) {
		_x = new UnlimitedFloatArray();
		_y = new UnlimitedFloatArray();
		_name = name;
		_xLabel = xlabel;
		_yLabel = ylabel;
	}
	public String getName() {
		return _name;
	}
	public String getTitle() {
		return _name;
	}
	public float getLastX() {
		return _x.getLastValue();
	}
	public String getXLabel() {
		return _xLabel;
	}
	public float getLastY() {
		return _y.getLastValue();
	}
	public String getYLabel() {
		return _yLabel;
	}
	/**
	 * @param x
	 *            The _x to set.
	 */
	public void setX(UnlimitedFloatArray x) {
		_x = x;
	}
	/**
	 * @param y
	 *            The _y to set.
	 */
	public void setY(UnlimitedFloatArray y) {
		_y = y;
	}
	/**
	 * @return Returns the _x.
	 */
	public UnlimitedFloatArray getXArray() {
		return _x;
	}
	/**
	 * @return Returns the _y.
	 */
	public UnlimitedFloatArray getYArray() {
		return _y;
	}
	/**
	 * @param _name
	 *            The _name to set.
	 */
	public void setName(String _name) {
		this._name = _name;
	}
	/**
	 * @param label
	 *            The _xLabel to set.
	 */
	public void setXLabel(String label) {
		_xLabel = label;
	}
	/**
	 * @param label
	 *            The _yLabel to set.
	 */
	public void setYLabel(String label) {
		_yLabel = label;
	}
}