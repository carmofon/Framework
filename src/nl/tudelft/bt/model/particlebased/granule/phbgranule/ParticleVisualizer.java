/*
 * Created on Sep 10, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Generation - Code and Comments
 */
package nl.tudelft.bt.model.particlebased.granule.phbgranule;

import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Collection;
import java.util.Iterator;

import javax.swing.*;

import org.jfree.chart.ChartPanel;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.util.Semaphore;

/**
 * @author jxavier
 * 
 * TODO To change the template for this generated type comment go to Window -
 * Preferences - Java - Code Generation - Code and Comments
 */
public class ParticleVisualizer extends ChartPanel {
	private Dimension _preferredSize = new Dimension(550, 300);
	private Color _bg = Color.WHITE;
	private float _x;
	private float _y;
	private float _w;
	private float _h;
	private float _hrect;
	private ParticleShape[] _bacteria;
	private Semaphore _sem = new Semaphore();
	private boolean _updated = true;
	private Ellipse2D.Float _circle = new Ellipse2D.Float();
	private ContourPlotHolder _cph;
	private Font _labelFont = new Font("Arial", Font.PLAIN, 12);

	class ParticleShape {
		private float _centerX;
		private float _centerY;
		private float _centerZ;
		private float _coreRadius;
		private Color _coreColor;

		/**
		 * create a new shape from the particle
		 * 
		 * @param b
		 */
		public ParticleShape(BiomassParticle b) {
			_centerX = b.getCenterX();
			_centerY = b.getCenterY();
			_centerZ = b.getCenterZ();
			_coreRadius = b.getCoreRadius();
			_coreColor = b.getColorCore();
		}

		/**
		 * Draw the shape
		 * 
		 * @param g
		 */
		public void draw(Graphics2D g) {
			float x = (_centerY - _coreRadius) / Model.model().systemSize.y
					* (float) _w + (float) _x;
			float y = -(_centerX + _coreRadius) / Model.model().systemSize.x
					* (float) _h + (float) _y;
			float w = (2 * _coreRadius) / Model.model().systemSize.y
					* (float) _w;
			float h = (2 * _coreRadius) / Model.model().systemSize.x
					* (float) _h;
			g.setPaint(_coreColor);
			_circle.setFrame(x, y, w, h);
			g.fill(_circle);
		}

	}

	/**
	 * Create the panel
	 */
	public ParticleVisualizer(ContourPlotHolder h) {
		// create the chart pannel with none of the fancy stuff
		super(h.getChart(), false, false, false, false, false);
		_cph = h;
		setMaximumDrawHeight(100000); //stop ChartPanel from scaling
		// output
		setMaximumDrawWidth(100000); //stop ChartPanel from scaling
		// output
		setHorizontalZoom(false);
		setVerticalZoom(false);
		setFillZoomRectangle(false);
		//Add a border
		setBorder(BorderFactory.createEtchedBorder());
		//
		setBackground(_bg);
		setOpaque(true);
	}

	public Dimension getPreferredSize() {
		return _preferredSize;
	}

	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		Graphics2D g2 = (Graphics2D) g;
		g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
				RenderingHints.VALUE_ANTIALIAS_ON);
		// get the available area for drawing
		Rectangle2D area = _cph.getDataArea();
		_w = (float) area.getWidth();
		_h = (float) area.getHeight();
		_x = (float) area.getX() + _w * 0.005f;
		_y = (float) area.getY() + _h;
		_hrect = 0.01f * _h;
		// draw the biofilm
		_sem.waitIfClosed();
		_sem.close();
		//don't draw the solid support
		//g2.setPaint(Color.BLACK);
		//g2.fill(new Rectangle2D.Float(_x, _y, _w, _hrect));
		//get the particles and draw
		if (_bacteria != null) {
			for (int i = 0; i < _bacteria.length; i++) {
				_bacteria[i].draw(g2);
			}
		}
		//set the font properties before rendering the text
		g2.setFont(_labelFont);
		g2.setPaint(Color.BLACK);
		//draw time label
		g2.drawString(formatTime(), _w * 0.05f, _h * 0.05f);
		//draw the scale
		g2.drawString("100 micron", _w * 0.80f, _h * 0.05f);		
		g2.fill(new Rectangle2D.Float( _w * 0.83f, _h * 0.06f, _w*0.1f, _hrect));
		_sem.open();
	}

	public void updateParticleData() {
		_sem.waitIfClosed();
		_sem.close();
		Collection bl = Model.model().getBiomassAsBiomassParticleCollection();
		int nbac = bl.size();
		_bacteria = new ParticleShape[nbac];
		int i = 0;
		for (Iterator iter = bl.iterator(); iter.hasNext();) {
			_bacteria[i++] = new ParticleShape((BiomassParticle) iter.next());
		}
		_cph.getOxygenConcentrations();
		_sem.open();
		repaint();
	}

	private String formatTime() {
		float t = Model.model().getTime();
		NumberFormat frm = new DecimalFormat("0.#");
		if (t > 24) {
			int day = (int) Math.floor(t / 24f);
			float h = t - ((float)day *24f);
			return new String("Time: " + day + " day " + frm.format(h) + " hour");			
		}
		return new String("Time: " + frm.format(t) + " h");
	}

}