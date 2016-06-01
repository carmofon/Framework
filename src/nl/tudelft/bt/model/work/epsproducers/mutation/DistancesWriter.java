package nl.tudelft.bt.model.work.epsproducers.mutation;

import java.util.Iterator;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.output.StateWriter;
import nl.tudelft.bt.model.exceptions.ModelException;
import nl.tudelft.bt.model.exceptions.ModelIOException;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

public class DistancesWriter extends StateWriter {
	static final private String DISTDIR = "distances";

	/*
	 * (non-Javadoc)
	 */
	public void write() throws ModelException {
		// write only if mutations already occured
		if (!MutatorBiomassParticleKeepDistances.mutationOccurred())
			return;
		try {
			String d = confirmSubDirectoryExists(DISTDIR);
			// file for distances of others
			String fn = d + "/distances"
					+ Model.model().getFormatedIterationNumber() + ".txt";
			java.io.File f = new java.io.File(fn);
			java.io.FileWriter fr = new java.io.FileWriter(f);
			// file for masses of the mutant
			String fn2 = d + "/mutant"
					+ Model.model().getFormatedIterationNumber() + ".txt";
			java.io.File f2 = new java.io.File(fn2);
			java.io.FileWriter fr2 = new java.io.FileWriter(f2);
			for (Iterator iter = Model.model().biomassContainer
					.getBiomassAsBiomassParticleCollection().iterator(); iter
					.hasNext();) {
				// get the current bacterium
				BiomassParticle b = (BiomassParticle) iter.next();
				if (b instanceof MutatorBiomassParticleKeepDistances) {
					MutatorBiomassParticleKeepDistances p = (MutatorBiomassParticleKeepDistances) b;
					if (!p.isMutant() & !p.isEpsOnly())
						fr.write("" + p.getDistanceToMutant() + "\t"
								+ p.getActiveMass() + "\n");
					else if (p.isMutant() & !p.isEpsOnly())
						fr2.write("" + p.getActiveMass() + "\n");
				}
			}
			fr.close();
			fr2.close();

		} catch (Exception e) {
			throw new ModelIOException("Error trying to write"
					+ " distances to file");
		}
	}
}
