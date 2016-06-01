package nl.tudelft.bt.model.povray;
import java.awt.Color;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
/**
 * Creates a union with a box (the biofilm carrier) and spheres (the bacteria)
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class Biofilm3D  implements Serializable {
	private Povray3DScene _pov;
	private VectorProperty translate;
	private VectorProperty rotate;
	private Box _biofilmCarrier;
	private ParticleWithCapsule[] _cells;
	private int _next;
	/**
	 * Construct a
	 */
	protected Biofilm3D(Povray3DScene pov) {
		_pov = pov;
		initializeAll();
	}
	protected Biofilm3D(Povray3DScene pov, int n) {
		_pov = pov;
		initializeAll();
		_cells = new ParticleWithCapsule[n];
	}
	private void initializeAll() {
		translate = new VectorProperty("translate");
		translate.setValues(-_pov.getX() * 0.4f, -_pov.getY() * 0.5f, -_pov
				.getZ() * 0.5f);
		rotate = new VectorProperty("rotate");
		rotate.setValues(0, 0, 90);
		_biofilmCarrier = new Box();
		_biofilmCarrier.setColor(0.2f, 0.2f, 0.2f);
		_biofilmCarrier.setCorner1(0f, 0f, 0f);
		_biofilmCarrier.setCorner2(-_pov.getY() * .01f, // thickness of box
				_pov.getY(), _pov.getZ());
	}
	/**
	 * Add a new cell
	 * 
	 * @param x
	 * @param y
	 * @param z
	 * @param rad
	 * @param r
	 * @param g
	 * @param b
	 */
	protected void addCell(float x, float y, float z, float rad, int r, int g,
			int b) {
		Color color = new Color(r, g, b);
		_cells[_next] = new ParticleWithCapsule();
		_cells[_next].setCenter(new ContinuousCoordinate(x, y, z));
		_cells[_next].setCoreRadius(rad);
		_cells[_next].setColorCore(color);
		_next++;
	}
	/**
	 * Writes the cell positions to pov file
	 * 
	 * @param f
	 * @throws IOException
	 */
	protected void toFile(FileWriter f) throws IOException {
		f.write("union {\n");
		f.write(_biofilmCarrier.toString());
		for (int i = 0; i < _cells.length; i++) {
			f.write(_cells[i].toString());
		}
		f.write("\t" + translate + "\n");
		f.write("\t" + rotate + "\n");
		f.write("}");
	}
	/**
	 * Writes the current cells position in model to file not using the _cells
	 * array, wich is too memory consuming
	 * 
	 * @param f
	 * @throws IOException
	 */
	protected void modelStateToFile(FileWriter f) throws IOException {
		biofilmHeaderToFile(f);
		particlesToFile(f);
		biofilmFooterToFile(f);
	}
	/**
	 * Writes the union open and the carrier to file
	 * 
	 * @param f
	 * @throws IOException
	 */
	protected void biofilmHeaderToFile(FileWriter f) throws IOException {
		f.write("union {\n");
		f.write(_biofilmCarrier.toString());
	}
	/**
	 * writes the union tranlations and rotations and close to file
	 * 
	 * @param f
	 * @throws IOException
	 */
	protected void biofilmFooterToFile(FileWriter f) throws IOException {
		f.write("\t" + translate + "\n");
		f.write("\t" + rotate + "\n");
		f.write("}");
	}
	/**
	 * Write the particles only to file
	 * 
	 * @param f
	 * @throws IOException
	 */
	protected void particlesToFile(FileWriter f) throws IOException {
		ArrayList bl = new ArrayList(Model.model()
				.getBiomassAsBiomassParticleCollection());
		for (Iterator iter = bl.iterator(); iter.hasNext();) {
			ParticleWithCapsule s = new ParticleWithCapsule(
					(BiomassParticle) iter.next());
			f.write(s.toString());
		}
	}
}
