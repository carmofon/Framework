/*
 * Created on 6-feb-2004 by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.detachment.levelset;

import java.io.Serializable;
import java.util.Collection;
import java.util.Iterator;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.util.UnlimitedFloatArray;

/**
 * A container that keeps trak of all biomass detached from the biomfilm.
 * Contains methods to return variable series of total detached biomass, biomass
 * detached by erosion mechanisms and biomass detached via sloughing.
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class DetachedBiomassContainer implements Serializable {
	private DetachedSeries _totalDetachedSeries;

	private DetachedSeries _erodedBiomassSeries;

	private DetachedSeries _sloughedBiomassSeries;

	/**
	 * 
	 */
	public DetachedBiomassContainer(UnlimitedFloatArray timeSeries) {
		_totalDetachedSeries = new DetachedSeries(this,
				"Total detached biomass", "Time [T]", "Detached biomass [M]",
				timeSeries);
		_erodedBiomassSeries = new DetachedSeries(this, "Total eroded biomass",
				"Time [T]", "Eroded biomass [M]", timeSeries);
		_sloughedBiomassSeries = new DetachedSeries(this,
				"Total sloughed biomass", "Time [T]", "Sloughed biomass [M]",
				timeSeries);
	}

	/**
	 * Add detached biomass to the eroded biomass and total detached biomass
	 * containers
	 * 
	 * @param c
	 */
	public void addToErodedBiomassSeries(BiomassSpecies.Composition c) {
		_totalDetachedSeries.add(c);
		_erodedBiomassSeries.add(c);
	}

	/**
	 * @param p
	 */
	public void addToDetachedBiomass(BiomassParticle p) {
		BiomassSpecies.Composition c = p.getComposition();
		_totalDetachedSeries.add(c);
		if (p.willDetachByErosion())
			_erodedBiomassSeries.add(c);
		else if (p.willDetachBySloughing())
			_sloughedBiomassSeries.add(c);
		// increment the calues of each particulate species in the composition
		Collection sp = Model.model().getParticulateSpecies();
		int speciesNumber = 0;
		for (Iterator iter = sp.iterator(); iter.hasNext();) {
			ParticulateSpecies element = (ParticulateSpecies) iter.next();
			element.addToDetachedBiomass(c.getSpeciesMass(element));
		}
	}

	/**
	 * @return Returns the time series of total detached biomass.
	 */
	public VariableSeries getTotalDetachedBiomassSeries() {
		return _totalDetachedSeries;
	}

	/**
	 * @return Returns the time series of eroded biomass.
	 */
	public DetachedSeries getErodedBiomassSeries() {
		return _erodedBiomassSeries;
	}

	/**
	 * @return Returns the time series of sloughed biomass.
	 */
	public DetachedSeries getSloughedBiomassSeries() {
		return _sloughedBiomassSeries;
	}
}