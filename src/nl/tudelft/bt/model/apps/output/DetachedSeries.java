/*
 * Created on Jun 30, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package nl.tudelft.bt.model.apps.output;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.detachment.levelset.DetachedBiomassContainer;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.util.UnlimitedFloatArray;

/**
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class DetachedSeries extends VariableSeries {
	private final DetachedBiomassContainer container;

	private ArrayList _detachedCompositions;

	/**
	 * @param name
	 *            name of series
	 * @param xlabel
	 * @param ylabel
	 */
	public DetachedSeries(DetachedBiomassContainer container, String name,
			String xlabel, String ylabel, UnlimitedFloatArray timeSeries) {
		super(name, xlabel, ylabel);
		this.container = container;
		setX(timeSeries);
		_detachedCompositions = new ArrayList();
	}

	/**
	 * Saves the composition and the time at which biomass detached
	 * 
	 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
	 */
	private class DetachedComposition {
		private float[] particulateMasses;

		private float time;

		/**
		 * @param c
		 *            composition of detached biomass
		 * @param t
		 *            time of detachment
		 */
		public DetachedComposition(BiomassSpecies.Composition c, float t) {
			time = t;
			addComposition(c);
		}

		/**
		 * Add a the masses of another compostion
		 * 
		 * @param c
		 */
		public void addComposition(BiomassSpecies.Composition c) {
			Collection sp = Model.model().getParticulateSpecies();
			// if the array is not initialized, initialize it
			if (particulateMasses == null)
				particulateMasses = new float[sp.size()];
			// add the species
			int speciesNumber = 0;
			for (Iterator iter = sp.iterator(); iter.hasNext();) {
				ParticulateSpecies element = (ParticulateSpecies) iter.next();
				particulateMasses[speciesNumber++] += c.getSpeciesMass(element);
			}
		}

		/**
		 * @return the sum of the masses of all species in this composition
		 */
		public float getTotalMass() {
			float m = 0;
			for (int i = 0; i < particulateMasses.length; i++) {
				m += particulateMasses[i];
			}
			return m;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.model.apps.output.VaribleSeries#getY()
	 */
	public float getLastY() {
		int sizeX = getXArray().getSize();
		int sizeY = getYArray().getSize();
		// if iterations are not synchronized with Y
		if (sizeY < sizeX) {
			// add all missing iterations
			for (int i = sizeY; i < sizeX; i++) {
				float t = getXArray().getValue(i);
				float m = 0;
				// iterate through the array, and remove any compositions
				// added to the series in the meanwhile
				int ncomp = _detachedCompositions.size();
				for (int j = 0; j < _detachedCompositions.size(); j++) {
					DetachedComposition d = (DetachedComposition) _detachedCompositions
							.get(j);
					m += d.getTotalMass();
					_detachedCompositions.remove(j);
					j--;
				}
				getYArray().add(m);
			}
		}
		// every time getY is invoked, the array is updated
		return super.getLastY();
	}

	/**
	 * Add an instance of composition corresponding to detached biomass
	 * 
	 * @param c
	 *            composition to add
	 */
	public void add(BiomassSpecies.Composition c) {
		// check if last detached composition added is at the same time
		int indexLast = _detachedCompositions.size() - 1;
		float time = Model.model().getTime();
		if (indexLast > 0) {
			DetachedComposition last = (DetachedComposition) _detachedCompositions
					.get(indexLast);
			if (time == last.time) {
				// if time is the same as as the composition added, only sum the
				last.addComposition(c);
				return;
			}

		}
		_detachedCompositions.add(new DetachedComposition(c, Model.model()
				.getTime()));
	}
}