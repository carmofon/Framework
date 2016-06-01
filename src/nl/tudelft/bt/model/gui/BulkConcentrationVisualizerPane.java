package nl.tudelft.bt.model.gui;

import java.util.Collection;
import java.util.Iterator;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.*;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;

import javax.swing.*;

/**
 * Applet to display biomass particles in present model state
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class BulkConcentrationVisualizerPane extends JPanel {
	// 
	private ApplicationComponent _application;
	private ChartFrame _frm;
	private XYSeriesCollection _data;
	private XYSeries[] _series1;
	private SoluteSpecies[] _soluteSpecies;
	/**
	 *  
	 */
	public BulkConcentrationVisualizerPane(ApplicationComponent c) {
		super();
		_application = c;
		// create some data...
		//create the data containers
		_data = new XYSeriesCollection();
		// create a chart...
		JFreeChart chart =
			ChartFactory.createXYLineChart(
				"Bulk concentration of solutes",
				"time [h]",
				"Concentration [g/l]",
				_data,
				PlotOrientation.VERTICAL,
				true,
				true,
				false);

		_frm = new ChartFrame("Bulk concentration of solutes", chart);
		_frm.pack();
		_frm.setVisible(true);
	}

	/**
	 * Update the displayed image showing present model state
	 */
	public void update() {
		// check if series are initialized, if not create series for each
		// colute species in the system
		if (_data.getSeriesCount() == 0)
			initializeSeries();
		else
			updateSerie();
		_frm.repaint();
	}

	/**
	 * Initialize the series for the plot by getting the solutespecies from
	 * model
	 */
	private void initializeSeries() {
		Collection s = Model.model().getSoluteSpecies();
		_soluteSpecies = new SoluteSpecies[s.size()];
		_series1 = new XYSeries[s.size()];
		int count = 0; // counter for array index
		for (Iterator iter = s.iterator(); iter.hasNext();) {
			SoluteSpecies element = (SoluteSpecies) iter.next();
			//copy colection to array
			_soluteSpecies[count] = element;
			//initialize a series for each solute species
			_series1[count] = new XYSeries(element.getName());
			// and add the first value
			float value = element.getBulkConcentration();
			_series1[count].add(Model.model().getTime(), value);
			_data.addSeries(_series1[count++]);
		}
	}

	/**
	 * Update the value of series with the value of the current iteration
	 */
	private void updateSerie() {
		for (int i = 0; i < _series1.length; i++) {
			_series1[i].add(
				Model.model().getTime(),
				_soluteSpecies[i].getBulkConcentration());
		}
	}
}
