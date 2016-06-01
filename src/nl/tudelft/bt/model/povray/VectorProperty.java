/*
 * Created on Jun 19, 2003
 * 
 * To change the template for this generated file go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
package nl.tudelft.bt.model.povray;

import java.io.Serializable;

/**
 * Vector property for PovRay descriptive language
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class VectorProperty  implements Serializable {
	private String _name;
	private float[] _values;

	public VectorProperty(String name) {
		_name = name;
	}

	public void setValues(float x, float y, float z) {
		_values = new float[3];
		_values[0] = x;
		_values[1] = y;
		_values[2] = z;
	}

	public void setValues(float x, float y, float z, float w) {
		_values = new float[4];
		_values[0] = x;
		_values[1] = y;
		_values[2] = z;
		_values[3] = w;
	}
	
	public String toString() {
		int n = _values.length;
		String out = _name + " <";
		for (int i = 0; i < n; i++) {
			out += " " + _values[i] + (i == n-1 ? " >" : ", ");
		}
		return out;
	}
}
