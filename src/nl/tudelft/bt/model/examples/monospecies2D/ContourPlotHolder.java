package nl.tudelft.bt.model.examples.monospecies2D;

import java.awt.Color;
import java.awt.GradientPaint;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.util.Date;


import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;

import org.jfree.chart.CrosshairInfo;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.ColorBar;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.ContourPlot;
import org.jfree.chart.plot.PlotRenderingInfo;
import org.jfree.chart.ui.ColorPalette;
import org.jfree.data.ContourDataset;
import org.jfree.data.DefaultContourDataset;
import org.jfree.ui.RectangleEdge;

/**
 *Class that creates the contour plot used by Monospecies2DGUI 
 *
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 * @see Monospecies2DGUI
 */
public class ContourPlotHolder {

	/** The x-axis. */
	private NumberAxis xAxis = null;

	/** The y-axis. */
	private NumberAxis yAxis = null;

	/** The z-axis. */
	private ColorBar zColorBar = null;

	/** The number of x values in the dataset. */
	private static int numX = Monospecies2DModel.GRIDSIZE;

	/** The number of y values in the dataset. */
	private static int numY = Monospecies2DModel.GRIDSIZE;

	private Double[] xValues;
	private Double[] yValues;
	private Double[] concentrationValues;

	/** The ratio. */
	private static double ratio = 0.0;

	/** The panel. */
	public ParticleVisualizer panel = null;
	private JFreeChart _chart;
	private DefaultContourDataset _data;

	private class ConcentrationPlot extends ContourPlot {
		private Rectangle2D _dataArea;
		/**
		 * @param dataset
		 * @param domainAxis
		 * @param rangeAxis
		 * @param colorBar
		 */
		public ConcentrationPlot(ContourDataset dataset, ValueAxis domainAxis,
				ValueAxis rangeAxis, ColorBar colorBar) {
			super(dataset, domainAxis, rangeAxis, colorBar);
		}
		public void render(Graphics2D g2, Rectangle2D dataArea,
				PlotRenderingInfo info, CrosshairInfo crosshairInfo) {
			super.render(g2, dataArea, info, crosshairInfo);
			_dataArea = dataArea;
		}
		/**
		 * @return Returns the _dataArea.
		 */
		public Rectangle2D getDataArea() {
			return _dataArea;
		}
	}
	/**
	 * Constructs a new demonstration application.
	 * 
	 * @param title
	 *            the frame title.
	 */
	public ContourPlotHolder() {
		_chart = createContourPlot();
		panel = new ParticleVisualizer(this);
	}

	/**
	 * Creates a ContourPlot chart.
	 * 
	 * @return the chart.
	 */
	private JFreeChart createContourPlot() {
		String title = "";
		String zAxisLabel = "Oxygen concentration [gO/l]";

		xAxis = new NumberAxis();
		xAxis.setAutoRangeIncludesZero(false);
		xAxis.setLowerMargin(0.0);
		xAxis.setUpperMargin(0.0);
		xAxis.setTickLabelsVisible(false);
		xAxis.setTickMarksVisible(false);

		yAxis = new NumberAxis();
		yAxis.setAutoRangeIncludesZero(false);
		yAxis.setLowerMargin(0.0);
		yAxis.setUpperMargin(0.0);
		yAxis.setTickLabelsVisible(false);
		yAxis.setTickMarksVisible(false);

		zColorBar = new ColorBar(zAxisLabel);
		ColorPalette pal = new BonePalette();
		zColorBar.setColorPalette(pal);
		zColorBar.getAxis().setTickMarksVisible(true);
				

		createDataset();
		ContourPlot plot = new ConcentrationPlot(_data, xAxis, yAxis, zColorBar);
		plot.setDataAreaRatio(ratio);
		plot.setColorBarLocation(RectangleEdge.RIGHT);

		JFreeChart chart = new JFreeChart(title, null, plot, false);
		chart.setBackgroundPaint(new GradientPaint(0, 0, Color.white, 0, 400,
				new Color(203, 197, 180)));

		return chart;

	}

	/**
	 * Creates a ContourDataset.
	 */
	private void createDataset() {
		int numValues = numX * numY;
		Date[] tmpDateX = new Date[numValues];
		double[] tmpDoubleX = new double[numValues];
		double[] tmpDoubleY = new double[numValues];

		xValues = new Double[numValues];
		yValues = new Double[numValues];
		concentrationValues = new Double[numValues];

		float stepX = Monospecies2DModel.SYSTEMSIZE / numX;
		float stepY = Monospecies2DModel.SYSTEMSIZE / numY;

		int c = 0;
		for (int i = 0; i < numX; i++) {
			for (int j = 0; j < numY; j++) {
				xValues[c] = new Double(i);
				yValues[c] = new Double(j);
				concentrationValues[c] = new Double(0);
				c++;
			}
		}
		getOxygenConcentrations();
		_data = null;

		_data = new DefaultContourDataset("Contouring", xValues, yValues,
				concentrationValues);

	}

	/**
	 * Gets the value of the oxygen concentrations from the model
	 */
	protected void getOxygenConcentrations() {
		try {
			//get the oxygen species
			SoluteSpecies oxygen = Model.model().getSoluteSpecies(0);
			float[][][] conc = oxygen.getFinestGrid();
			int c = 0;
			for (int i = 0; i < numX; i++) {
				for (int j = 0; j < numY; j++) {
					concentrationValues[c] = new Double(conc[j+1][i+1][1]);
					c++;
				}
			}
			_data.initialize(xValues, yValues, concentrationValues);
			((ContourPlot) _chart.getPlot()).setDataset(_data);
		} catch (NullPointerException e) {
			// Happens if system is not yet initialized
			// do nothing, all values are set to 0
		}
	}
	/**
	 * @return Returns the _chart.
	 */
	public JFreeChart getChart() {
		return _chart;
	}
	
	/**
	 * Get the area used for ploting the data
	 * 
	 * @return rectangle with area used for ploting data
	 */
	public Rectangle2D getDataArea() {
		return ((ConcentrationPlot)_chart.getPlot()).getDataArea();
	}
}