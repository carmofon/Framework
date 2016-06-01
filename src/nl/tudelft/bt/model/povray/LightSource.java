package nl.tudelft.bt.model.povray;

import java.io.Serializable;

/**
 * Light source for povray scene
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class LightSource  implements Serializable {
	private VectorProperty location;
	private VectorProperty color;

	/**
	 * 
	 */
	public LightSource() {
		location = new VectorProperty("");
		color = new VectorProperty("color rgb");
	}

	/**
	 * @param fs
	 */
	public void setLocation(float x, float y, float z) {
		location.setValues(x, y, z);
	}

	/**
	 * @param fs
	 */
	public void setColor(float r, float g, float b) {
		color.setValues(r, g, b);
	}

	public String toString() {
		return "light_source {\n"
			+ "\t "
			+ location
			+ "\n\t"
			+ color
			+ "\n"
			+ "}\n";
	}
}
