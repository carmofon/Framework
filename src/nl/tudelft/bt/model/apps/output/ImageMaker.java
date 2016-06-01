package nl.tudelft.bt.model.apps.output;

import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.util.ColorMaps;

import java.awt.*;
import java.util.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.DataBuffer;
import java.awt.image.DataBufferInt;
import java.awt.image.Raster;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.*;

import ch.randelshofer.media.quicktime.QuickTimeWriter;

/**
 * Applet to display biomass particles in present model state
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ImageMaker {
	//
	static final int WIDTH = 500;
	static final int HEIGHT = 500;

	final static Color bg = Color.white;
	final static Color fg = Color.black;

	final static BasicStroke stroke = new BasicStroke(1.0f);

	private ParticleShape[] _bacteria;
	private float _height;
	private float _width;
	private float _frac; // deepness ratio
	private Dimension _d;

	private QuickTimeWriter _quickTimeWriter;
	private static boolean _writing = false;

	private SoluteSpecies _solute;

	private BufferedImage bi;
	private Graphics2D ig2;

	class ParticleShape {
		private ContinuousCoordinate _center;
		private float _capsuleRadius;
		private float _coreRadius;
		private Color _capsuleColor;
		private Color _coreColor;

		/**
		 * create a new shape from the particle
		 * 
		 * @param b
		 */
		public ParticleShape(BiomassParticle b) {
			_center = b.getCenter();
			_capsuleRadius = b.getRadius();
			_coreRadius = b.getCoreRadius();
			_capsuleColor = b.getColorCapsule();
			_coreColor = b.getColorCore();
		}

		/**
		 * Draw the shape
		 * 
		 * @param g2
		 */
		public void draw(Graphics2D g2) {
			g2.setPaint(_capsuleColor);
			g2.fill(bacteriaCapsuleCircle(_center, _capsuleRadius));
			g2.setPaint(_coreColor);
			g2.fill(bacteriaCoreCircle(_center, _coreRadius));
		}
	}

	/**
	 *  
	 */
	public ImageMaker(SoluteSpecies s, String movieFileName) throws IOException {
		_d = new Dimension(WIDTH, HEIGHT);
		_quickTimeWriter = new QuickTimeWriter(new File(movieFileName));
		_quickTimeWriter.addVideoTrack(QuickTimeWriter.VideoFormat.PNG, 30,
				WIDTH, HEIGHT);
		_solute = s;

		// start thread to catch program shutdown
		Runtime.getRuntime().addShutdownHook(new Thread() {
			@Override
			public void run() {
				try {
					while (_writing) {
						Thread.sleep(100);
						System.out
								.println("couldn' write movie, trying again...");
					}
					_quickTimeWriter.close();
				} catch (Exception e) {
					e.printStackTrace();
				}
				System.out.println("Simulation done! Closing quicktime movie.");
			}
		});
	}

	/**
	 * Update the displayed image showing present model state
	 */
	public void update() {
		ArrayList bl = new ArrayList(Model.model()
				.getBiomassAsBiomassParticleCollection());
		Collections.sort(bl, new BiomassParticle.DepthComparator());
		int nbac = bl.size();
		_bacteria = new ParticleShape[nbac];
		for (int i = 0; i < nbac; i++) {
			BiomassParticle b = (BiomassParticle) bl.get(i);
			_bacteria[i] = new ParticleShape(b);
		}
	}

	/**
	 * Return a circular shape with the position and radius of particle to draw
	 * 
	 * @param b
	 *            particle to draw
	 * @return corresponding shape to be drawn
	 */
	private Shape bacteriaCoreCircle(ContinuousCoordinate c, float r) {
		float x = (c.y - c.z * 0.5f - r) / _width;
		x = (float) ((0.1 + _frac + (0.8 - _frac) * x) * _d.width);
		float y = (_height - c.x - c.z * 0.5f - r) / _height;
		y = (float) ((0.1 + _frac + (0.8 - _frac) * y) * _d.height);
		float w = (2 * r) / _width;
		w = (float) (0.8 * w * _d.width);
		float h = (2 * r) / _height;
		h = (float) (0.8 * h * _d.height);
		return new Ellipse2D.Float(x, y, w, h);
	}

	private Shape bacteriaCapsuleCircle(ContinuousCoordinate c, float r) {
		float x = (c.y - c.z * 0.5f - r) / _width;
		x = (float) ((0.1 + _frac + (0.8 - _frac) * x) * _d.width);
		float y = (_height - c.x - c.z * 0.5f - r) / _height;
		y = (float) ((0.1 + _frac + (0.8 - _frac) * y) * _d.height);
		float w = (2 * r) / _width;
		w = (float) (0.8 * w * _d.width);
		float h = (2 * r) / _height;
		h = (float) (0.8 * h * _d.height);
		return new Ellipse2D.Float(x, y, w, h);
	}

	/**
	 * Draw representation of solid substratum. Diferent implementation if
	 * system is 2D or 3D
	 * 
	 * @param g
	 */
	private void drawSolidSurface(Graphics2D g) {
		float x;
		float y;
		float w;
		float h;
		// draw rectangle representing solid substratum
		g.setPaint(fg);
		g.setStroke(stroke);
		if (Model.model().getDimensionality() == 2) {
			// 2D case
			x = 0.10f * _d.width;
			y = 0.90f * _d.height;
			w = 0.8f * _d.width;
			h = 0.05f * _d.height;
			g.draw(new Rectangle2D.Float(x, y, w, h));
		} else {
			// 3D case
			x = (0.10f + _frac) * _d.width;
			y = 0.90f * _d.height;
			w = (0.8f - _frac) * _d.width;
			h = 0.05f * _d.height;
			// draw trapezoid
			float x1Points[] = { x, x + w, x + _d.width * _frac,
					x - _d.width * _frac };
			float y1Points[] = { y, y, y - _d.height * _frac,
					y - _d.height * _frac };
			GeneralPath polygon = new GeneralPath(GeneralPath.WIND_EVEN_ODD,
					x1Points.length);
			polygon.moveTo(x1Points[0], y1Points[0]);

			for (int index = 1; index < x1Points.length; index++) {
				polygon.lineTo(x1Points[index], y1Points[index]);
			}
			;

			polygon.closePath();
			g.draw(polygon);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.awt.Component#paint(java.awt.Graphics)
	 */
	public void writeFrame(String fn) throws IOException {
		setSystemSize();
		update();

		// TYPE_INT_ARGB specifies the image format: 8-bit RGBA packed
		// into integer pixels
		bi = new BufferedImage(_d.width, _d.height, BufferedImage.TYPE_INT_ARGB);

		ig2 = bi.createGraphics();
		ig2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
				RenderingHints.VALUE_ANTIALIAS_ON);

		// clear previous image
		ig2.setBackground(Color.white);
		ig2.clearRect(0, 0, _d.width, _d.height);

		// write the gradients
		if ((Model.model().getDimensionality() == 2) && (_solute != null)) {
			float[][][] soluteMatrix = _solute.getFinestGrid();
			float ratioWidth = ((float) _d.width)
					/ ((float) (soluteMatrix[0].length - 2));
			float ratioHeight = ((float) _d.height)
					/ ((float) (soluteMatrix.length - 2));
			for (int i = 1; i < soluteMatrix.length - 2; i++) {
				for (int j = 1; j < soluteMatrix[0].length - 2; j++) {
					float x = ((float) j - 1) * ratioWidth * 0.8f + 0.1f
							* _d.width;
					float y = ((float) (soluteMatrix.length - i - 2))
							* ratioHeight * 0.8f + 0.1f * _d.height;
					float w = ratioWidth;
					float h = ratioHeight;
					float max = _solute.getMaximumValue();
					float min = _solute.getMinimumValue();
					float v = (soluteMatrix[i][j][1] - min) / (max-min);
					ig2.setPaint(ColorMaps.getBoneColor(v));
					ig2.fill(new Rectangle2D.Float(x, y, w, h));
				}
			}
		}
		// draw the biomass
		try {

			// draw solid substratum
			// drawSolidSurface(ig2);

			// draw particles
			if (_bacteria != null)
				for (int i = 0; i < _bacteria.length; i++) {
					_bacteria[i].draw(ig2);
				}
			// ig2.dispose();
			ImageIO.write(bi, "PNG", new File(fn + ".png"));
			_writing = true;
			_quickTimeWriter.writeFrame(0, bi, 2);
			_writing = false;
		} catch (NullPointerException e) {
			e.printStackTrace();
			// do nothing
			// the purpose of this catch is for when visualizer is
			// initialized before the system size
			return;
		}
		ig2.dispose();
	}

	/**
	 * (Re)define the system size atrributes
	 */
	public void setSystemSize() {
		// set the system size:
		_height = Model.model().systemSize.y;
		_width = Model.model().systemSize.x;
		// for 2D perspective projection
		if (Model.model().systemSize.z > 0)
			// if system is 3D
			_frac = 0.8f / (1.0f + Model.model().systemSize.x
					/ Model.model().systemSize.z * 2.0f);
		else
			// if system is 2D
			_frac = 0;
	}

}
