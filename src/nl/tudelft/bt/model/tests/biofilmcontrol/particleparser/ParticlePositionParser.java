package nl.tudelft.bt.model.tests.biofilmcontrol.particleparser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.ModelHandler;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * Implementes a class that may be used to read particle postition files and
 * used these positions to inoculate the system. This way a simulation may be
 * restarted from a given pre-computed state
 * 
 * @author jxavier
 */
public class ParticlePositionParser {
	private ArrayList _particleList;

	private BiomassSpecies _biomassSpecies;

	private float _coreSpecificMass;

	private ParticulateSpecies _coreSpecies;

	private float _capsuleSpecificMass;

	private ParticulateSpecies _capsuleSpecies;

	private class ParticleInformation {
		private float x;

		private float y;

		private float z;

		private float[] masses;

		/**
		 * Create an instance by parsing the tokens contained in a line of the
		 * particle position file
		 * 
		 * @param tokens
		 */
		private ParticleInformation(String[] tokens) {
			x = Float.parseFloat(tokens[0]);
			y = Float.parseFloat(tokens[1]);
			z = Float.parseFloat(tokens[2]);
			// compute masses of core and particle from the volume of core and
			// capsule
			float radiusCore = Float.parseFloat(tokens[3]);
			float radiusCapsule = Float.parseFloat(tokens[7]);
			float volumeCore = (Model.model().getDimensionality() == 3
					? ExtraMath.volumeOfASphere(radiusCore)
					: ExtraMath.volumeOfACylinder(radiusCore, Model.model()
							.get2DSystem3rdDimension()));
			float volumeCapsule = (Model.model().getDimensionality() == 3
					? ExtraMath.volumeOfASphere(radiusCapsule)
					: ExtraMath.volumeOfACylinder(radiusCapsule, Model.model()
							.get2DSystem3rdDimension()))
					- volumeCore;
			float massCore = volumeCore * _coreSpecificMass;
			float massCapsule = volumeCapsule * _capsuleSpecificMass;
			masses = new float[2];
			masses[0] = massCore;
			masses[1] = massCapsule;
		}
	}

	/**
	 * @param fileName
	 * @param sp
	 * @param coreSpecificMass
	 * @param coreSpecies
	 *            The species constituting the core
	 * @param capsuleSpecificMass
	 * @param capsuleSpecies
	 *            The species contituting the capsule
	 */
	public ParticlePositionParser(String fileName, BiomassSpecies sp,
			float coreSpecificMass, ParticulateSpecies coreSpecies,
			float capsuleSpecificMass,
			nl.tudelft.bt.model.multigrid.ParticulateSpecies capsuleSpecies) {
		// set attributes
		_biomassSpecies = sp;
		_coreSpecificMass = coreSpecificMass;
		_coreSpecies = coreSpecies;
		_capsuleSpecificMass = capsuleSpecificMass;
		_capsuleSpecies = capsuleSpecies;
		// read and parse the file
		File d = new File(fileName);
		try {
			BufferedReader br = new BufferedReader(new FileReader(d));
			String line;
			String[] tokens;
			_particleList = new ArrayList();
			while ((line = br.readLine()) != null) {
				tokens = line.split("\t");
				ParticleInformation pi = new ParticleInformation(tokens);
				_particleList.add(pi);
			}
		} catch (FileNotFoundException e) {
			throw new ModelRuntimeException("file " + fileName + " not found.");
		} catch (IOException e) {
			throw new ModelRuntimeException(e.toString());
		}
	}

	/**
	 * Initiale the system state by palcing biomass at each of the locations
	 * read from the file
	 * 
	 * @param h
	 */
	public void initializeSystemState(ModelHandler h) {
		for (Iterator iter = _particleList.iterator(); iter.hasNext();) {
			ParticleInformation pi = (ParticleInformation) iter.next();
			BiomassParticle b = h.placeBiomass(_biomassSpecies, pi.x, pi.y,
					pi.z);
			b.setMass(_coreSpecies, pi.masses[0]);
			b.setMass(_capsuleSpecies, pi.masses[1]);
		}
	}
}