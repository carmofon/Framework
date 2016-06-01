package nl.tudelft.bt.model.gui;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.output.VariableSeries;
import javax.swing.*;
/**
 * Applet to display biomass particles in present model state
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class VariableSeriesVisualizerPane extends JPanel {
	// 
	private ApplicationComponent _application;
	private ChartFrame _frm;
	private XYSeriesCollection _data;
	private XYSeries[] _series;
	private VariableSeries[] _var;
	/**
	 *  
	 */
	public VariableSeriesVisualizerPane(ApplicationComponent c,
			VariableSeries[] v, String title) {
		super();
		_application = c;
		_var = v;
		// create some data...
		//create the data containers
		_data = new XYSeriesCollection();
		// create a chart...
		JFreeChart chart = ChartFactory.createXYLineChart(title,
				_var[0].getXLabel(), _var[0].getYLabel(), _data,
				PlotOrientation.VERTICAL, true, true, false);
		_frm = new ChartFrame(title, chart);
		_frm.pack();
		_frm.setVisible(true);
		//
		_series = new XYSeries[_var.length];
	}
	/**
	 * Update the displayed image showing present model state
	 */
	public void update() {
		// check if series are initialized, if not create series for each
		// colute species in the system
		if (_data.getSeriesCount() == 0)
			initializeSeries();
		for (int i = 0; i < _var.length; i++)
			//update variables if not updated yet
			if (_series[i].getItemCount() < _var[i].getXArray().getSize())
				_series[i].add(_var[i].getLastX(), _var[i].getLastY());
		_frm.repaint();
	}
	/**
	 * Initialize the series for the plot by getting the solutespecies from
	 * model
	 */
	private void initializeSeries() {
		for (int i = 0; i < _var.length; i++) {
			_series[i] = new XYSeries(_var[i].getName());
			_data.addSeries(_series[i]);
		}
	}
}
