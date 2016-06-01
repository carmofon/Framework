package nl.tudelft.bt.model.povray;

import java.io.Serializable;

/**
 * Box object
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class Box  implements Serializable {
	private VectorProperty corner1;
	private VectorProperty corner2;
	private VectorProperty color;

	public Box() {
		corner1 = new VectorProperty("");
		corner2 = new VectorProperty("");
		color = new VectorProperty("color rgb");
	}

	/**
	 * @param fs
	 */
	public void setColor(float r, float g, float b) {
		color.setValues(r, g, b);
	}

	/**
	 * @param fs
	 */
	public void setCorner1(float x, float y, float z) {
		corner1.setValues(x, y, z);
	}

	/**
	 * @param fs
	 */
	public void setCorner2(float x, float y, float z) {
		corner2.setValues(x, y, z);
	}

	public String toString() {
		return "box {\n"
			+ "\t "
			+ corner1
			+ "\n"
			+ "\t "
			+ corner2
			+ "\n"
			+ "\t pigment { "
			+ color
			+ " }\n"
			+ "\t\tfinish {\n"
			+ "\t\t\t phong 0.9\n"
			+ "\t\t\t phong_size 60\n"
			+ "\t\t metallic }\n"
			+ "}\n";
	}
}
