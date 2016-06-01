package nl.tudelft.bt.model.povray;

import java.io.Serializable;

/**
 * The scene background
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class Background  implements Serializable {
	private VectorProperty color;

	public Background() {
		color = new VectorProperty("color rgb");
	}

	/**
	 * @param fs
	 */
	public void setColor(float r, float g, float b) {
		color.setValues(r, g, b);
	}

	public String toString() {
		return "background {\n"
			+ "\t" + color + "\n"
			+ "}\n";
	}
}
