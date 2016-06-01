/*
 * File created originally on Mar 12, 2006
 */
package nl.tudelft.bt.model.particlebased.tube.detachment;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.DiscreteCoordinate;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.detachment.DetachmentHandler;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.multigrid.BiomassBooleanMatrix;
import nl.tudelft.bt.model.multigrid.MultigridUtils;
import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.multigrid.TubeVelocityField;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.particlebased.tube.TubeBiomassParticleContainer;
import nl.tudelft.bt.model.util.ExtraMath;
import nl.tudelft.bt.model.work.tube.MonoSpeciesTube;

/**
 * Includes the algorithms for determining which pices of biofilm will detach as
 * a consequence of the shear stress in the tube
 * 
 * @author jxavier
 */
public class DetachmentInducedByShear implements DetachmentHandler {
	private ArrayList _border;

	private SetOfGridNodes _borderWithValue;

	private Segments _segmentsToDetach;

	private BiomassBooleanMatrix _biomass;

	private TubeBiomassParticleContainer _tube;

	private float _breakingShearStress;

	private static final float NSTEPSR = 50f;

	private boolean[][][] _matrixOfLocationsToDetach;

	/**
	 * Part of the testing phase for tube reactors, implements a grid node in 2D
	 * coordinates (Tube2D)
	 * 
	 * @author jxavier
	 */
	private class GridNode2DWithValue extends DiscreteCoordinate {
		float value;

		float x;

		float y;

		float r;

		float theta;

		float xCentered;

		float yCentered;

		/**
		 * Build new instance with value set to 0
		 * 
		 * @param i
		 *            i coordinate
		 * @param j
		 *            j coordinate
		 */
		GridNode2DWithValue(int i, int j) {
			this.i = i;
			this.j = j;
			k = 1;
			// convert matrix indexes to continuous coordinates in cartesian
			// and polar coordinate system
			ContinuousCoordinate c = MultigridVariable
					.convertDiscreteToContinuous(this);
			x = c.x;
			y = c.y;
			r = _tube.rConvertToPolar(x, y);
			theta = _tube.thetaConvertToPolar(x, y);
			xCentered = x - _tube.getHalfSystemSide();
			yCentered = y - _tube.getHalfSystemSide();
			value = 0; // default value
		}

		/**
		 * Initialize with k = 2;
		 * 
		 * @param c
		 *            the continuous coordinate to copy here
		 * @param v
		 *            the value (e.g. the shear stress) associated with location
		 */
		GridNode2DWithValue(DiscreteCoordinate d, float v) {
			this(d.i, d.j);
			value = v;
		}

		/**
		 * @param n
		 * @return distance to point defined by n
		 */
		float distance(GridNode2DWithValue n) {
			return ExtraMath.pointDistance(x, y, n.x, n.y);
		}

		/**
		 * @param n
		 * @return the average value between this point and n
		 */
		float averageShear(GridNode2DWithValue n) {
			return (value + n.value) * 0.5f;
		}

		/**
		 * @param n
		 * @return distance to point defined by n
		 */
		float rotation(GridNode2DWithValue n) {
			return (theta - n.theta + (float) (2 * Math.PI))
					% ((float) (2 * Math.PI));
		}

		/**
		 * @param i
		 * @param j
		 * @return rue if location is the same
		 */
		boolean isSameLocation(int i, int j) {
			return ((this.i == i) & (this.j == j));
		}

		public boolean equals(Object obj) {
			return isSameLocation(((GridNode2DWithValue) obj).i,
					((GridNode2DWithValue) obj).j);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Object#toString()
		 */
		public String toString() {
			return "" + i + "\t" + j + "\t" + x + "\t" + y + "\t" + xCentered
					+ "\t" + yCentered + "\t" + r + "\t" + theta + "\t" + value
					+ "\n";
		}
	}

	/**
	 * Implements a list of grid nodes
	 * 
	 * @author jxavier
	 */
	private class SetOfGridNodes {
		ArrayList _list = new ArrayList();

		float archDistance = 0;

		float pathDistance = 0;

		float ratio;

		float shear;

		float radius;

		/**
		 * default constructor
		 */
		public SetOfGridNodes() {
		}

		/**
		 * Copy constructor (and add radius)
		 */
		public SetOfGridNodes(SetOfGridNodes s, float r) {
			for (Iterator iter = s._list.iterator(); iter.hasNext();) {
				GridNode2DWithValue element = (GridNode2DWithValue) iter.next();
				_list.add(element);
			}
			radius = r;
		}

		/**
		 * Constructor that only initializes the value of the radius
		 * 
		 * @param r
		 */
		public SetOfGridNodes(float r) {
			radius = r;
		}

		/**
		 * set the value of the radius for this segment
		 * 
		 * @param r
		 */
		void setRadius(float r) {
			radius = r;
		}

		/**
		 * add a node to the segment
		 * 
		 * @param n
		 */
		void add(GridNode2DWithValue n) {
			_list.add(n);
		}

		/**
		 * Clear the contents of the _list and set values of attributes to 0
		 */
		void clearContents() {
			_list.clear();
			archDistance = 0;
			pathDistance = 0;
			ratio = 0;
			shear = 0;
			radius = 0;
		}

		/**
		 * Pass the contents to n and clear the contents
		 * 
		 * @param n
		 */
		void passContentsTo(SetOfGridNodes n) {
			n.addContentsOf(this);
			this.clearContents();
		}

		/**
		 * add the nodes in a set n to this segment.
		 * 
		 * @param n
		 */
		void addContentsOf(SetOfGridNodes n) {
			for (int i = 0; i < n.size(); i++) {
				_list.add(n.get(i));
			}
		}

		/**
		 * get node at location i
		 * 
		 * @param i
		 * @return node at location i
		 */
		GridNode2DWithValue get(int i) {
			return (GridNode2DWithValue) _list.get(i);
		}

		/**
		 * @return the number of nodes in the segment
		 */
		int size() {
			return _list.size();
		}

		/**
		 * Compute the average shear applied on the arch
		 */
		void computeShear() {
			// compute the path distance and the path integral of shear
			pathDistance = 0;
			shear = 0;
			for (int i = 0; i < size() - 1; i++) {
				float d = get(i).distance(get(i + 1));
				pathDistance += d;
				shear += get(i).averageShear(get(i + 1)) * d;
			}
			// compute the arch distance
			if (radius == 0)
				throw new ModelRuntimeException(
						"Trying to compute arch for segment"
								+ " with no radius value");
			// testing code
			GridNode2DWithValue p1 = get(0);
			GridNode2DWithValue p2 = get(size() - 1);
			// determine the angle
			float angle = p1.rotation(p2);
			// compute the arch distance between the two points
			archDistance = angle * radius;
			// finish computing the average, dividing by the path length
			shear /= archDistance;
		}

		public String toString() {
			String out = "";
			// header - uncomment following code
			// out += "Segment properties:\n";
			// out += "path distance: " + pathDistance + "; arch distance:"
			// + archDistance + "; shear: " + shear + "\n";
			// out += "Nodes:\n";
			for (Iterator iter = _list.iterator(); iter.hasNext();) {
				GridNode2DWithValue element = (GridNode2DWithValue) iter.next();
				out += element.toString();
			}
			return out;
		}

		/**
		 * @param i
		 * @param j
		 * @return true if location is contained, false otherwise
		 */
		public boolean contains(int i, int j) {
			for (int n = 0; n < this.size(); n++) {
				if (this.get(n).isSameLocation(i, j))
					return true;
			}
			return false;
		}

		/**
		 * Check if node is contained in the set
		 * 
		 * @param n
		 * @return true is node is contained, false otherwise
		 */
		public boolean contains(GridNode2DWithValue n) {
			return contains(n.i, n.j);
		}

		/**
		 * @param smaller
		 * @return true if all elements of smaller are contained in this segment
		 */
		public boolean contains(SetOfGridNodes smaller) {
			if (smaller.size() > size())
				// don't even waste time checking
				return false;
			// iterate through smaller
			for (int i = 0; i < smaller.size(); i++)
				if (!contains(smaller.get(i)))
					// found a point in smaller which is not in this
					return false;
			// cheked all elements in smaller are found here
			return true;
		}
	}

	/**
	 * The segment to detach will also store the locations in the grid that will
	 * be detaching
	 * 
	 * @author jxavier
	 */
	private class SegmentToDetach extends SetOfGridNodes {
		SetOfGridNodes _locationsToDetach;

		/**
		 * Initialize a segment to detach with the data from a set that was
		 * evaluated
		 * 
		 * @param s
		 */
		public SegmentToDetach(SetOfGridNodes s) {
			super(s, s.radius);
			// find the locations to detach
			_locationsToDetach = new SetOfGridNodes(radius);
			SetOfGridNodes trial = new SetOfGridNodes(radius);
			// initialize found with s so that the first test is passed
			SetOfGridNodes found = s;
			// start detection
			while (found.size() > 0) {
				// pass the contents of found to trial
				found.passContentsTo(trial);
				// detect
				for (int n = 0; n < trial.size(); n++) {
					for (int i = -1; i <= 1; i++) {
						for (int j = -1; j <= 1; j++) {
							int iPoint = trial.get(n).i + i;
							int jPoint = trial.get(n).j + j;
							// add the point if it is a non neighbor and if it
							// was not found before
							if (_biomass.isCarrierOrBiomass2D(iPoint, jPoint)
									& (!_locationsToDetach.contains(iPoint,
											jPoint))
									& (!trial.contains(iPoint, jPoint))
									& (!found.contains(iPoint, jPoint))
									& isInsideRadius(radius, iPoint, jPoint)) {
								found.add(new GridNode2DWithValue(iPoint,
										jPoint));
							}
						}

					}
				}
			}
			trial.passContentsTo(_locationsToDetach);
		}

		/**
		 * @param r
		 * @param i
		 * @param j
		 * @return true if the radius corresponding to (i,j) is greater or equal
		 *         to r
		 */
		private boolean isInsideRadius(float r, int i, int j) {
			float x = MultigridVariable.convertDiscreteToContinuousX(i);
			float y = MultigridVariable.convertDiscreteToContinuousX(j);
			return (_tube.rConvertToPolar(x, y) <= r);
		}

		/**
		 * Call setLocationsToDetach for each of the segments in the this set of
		 * segments
		 * 
		 * @param matrix
		 */
		void setLocationsToDetach(boolean[][][] matrix) {
			for (int i = 0; i < _locationsToDetach.size(); i++) {
				GridNode2DWithValue node = ((GridNode2DWithValue) _locationsToDetach
						.get(i));
				matrix[node.i][node.j][1] = true;
			}
		}
	}

	/**
	 * Implements a list of Segment
	 * 
	 * @author jxavier
	 */
	private class Segments {
		ArrayList _list = new ArrayList();

		/**
		 * Add if shear of segment is above the breaking shear
		 * 
		 * @param s
		 */
		void addIfAboveTheBreakingShear(SetOfGridNodes s) {
			s.computeShear();
			// check if segment should detach
			if (s.shear > _breakingShearStress)
				addIfNotContainedInOtherSegments(s);
		}

		/**
		 * iterate through the list and check is any segment in the list already
		 * contains the segment
		 * 
		 * @param s
		 *            only add segments that
		 */
		private void addIfNotContainedInOtherSegments(SetOfGridNodes s) {
			for (int i = 0; i < size(); i++) {
				if (get(i).contains(s))
					return; // return without adding to list
			}
			// if so, convert it to a SegmentToDetach (which also
			// computes which grid nodes inside the segment will detach
			_list.add(new SegmentToDetach(s));
		}

		/**
		 * @param i
		 * @return segment i
		 */
		SegmentToDetach get(int i) {
			return (SegmentToDetach) _list.get(i);
		}

		/**
		 * @return size of the array
		 */
		int size() {
			return _list.size();
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Object#toString()
		 */
		public String toString() {
			String out = "";
			for (Iterator iter = _list.iterator(); iter.hasNext();) {
				SetOfGridNodes s = (SetOfGridNodes) iter.next();
				out += "Segment begin ------------:\n";
				out += s.toString();
				out += "end ------------:\n";
			}
			return out;
		}

		/**
		 * Detect which are the segments to detach from _borderWithValue
		 * 
		 * @param r
		 *            the radius at which segments to detach are being computed
		 */
		void detectSegmentsToDetach(float r) {
			int nPoints = _borderWithValue.size();
			// find the first grid outside radius
			int i = 0;
			for (; i < nPoints + 1; i++) {
				if (i == nPoints) {
					// for the case that all points are inside the radius
					// the segments are all in one single segment
					this.addIfAboveTheBreakingShear(new SetOfGridNodes(
							_borderWithValue, r));
					return;
				}
				// if pioint is outside radius then break
				if (_borderWithValue.get(i).r > r)
					break;
			}
			// keep the i value, since it will be used to detect if algorithm
			// was
			// completed
			int iStart = i;
			// start detecting segments
			while (true) {
				// travel allong path until first point inside r is found
				while (_borderWithValue.get(i).r > r) {
					i = (i + 1) % nPoints;
					// check if cycle is terminated
					if (i == iStart)
						return;
				}
				// colect all points until another point outside r is found
				SetOfGridNodes segment = new SetOfGridNodes(r);
				while (_borderWithValue.get(i).r <= r) {
					segment.add(_borderWithValue.get(i));
					i = (i + 1) % nPoints;
				}
				// add only if size is greater that 1
				if (segment.size() > 1) {
					this.addIfAboveTheBreakingShear(segment);
				}
				// check if cycle is terminated
				if (i == iStart)
					return;
			}
		}

		/**
		 * Call setLocationsToDetach for each of the segments in the this set of
		 * segments
		 * 
		 * @param matrix
		 */
		void setLocationsToDetach(boolean[][][] matrix) {
			for (int i = 0; i < size(); i++) {
				((SegmentToDetach) (get(i))).setLocationsToDetach(matrix);
			}
		}
	}

	/**
	 * Initialize the detachment handler
	 */
	public DetachmentInducedByShear(float detachmentShearStress) {
		_tube = (TubeBiomassParticleContainer) (Model.model().biomassContainer);
		_breakingShearStress = detachmentShearStress;
		_matrixOfLocationsToDetach = MultigridVariable
				.create3DBooleanMatrixWithFinnerResolution();
	}

	public void detach() {
		// get the updated matrix
		_biomass = Model.model().getBiomassBooleanMatrix();
		_border = _biomass.findTheBorderPointsInTube();
		// Compute the shear in each point;
		_borderWithValue = computeTheAverageShearForEachPointInBorder(_border,
				_biomass);
		// for a a range of radius, find all the segments
		// that have a shear applied which is greater than the breaking
		// shear
		float rTube = _tube.getRadiusOfTubularReactor() * 0.9f;
		float rStep = rTube * 1.1f / NSTEPSR;
		// reset the segments to detach
		_segmentsToDetach = new Segments();
		for (float r = rTube * 1.5f; r > 0; r -= rStep) {
			// add progressivelly to the list of segments to detach
			_segmentsToDetach.detectSegmentsToDetach(r);
		}
		System.out
				.println("Detached " + _segmentsToDetach.size() + " segments"); // TODO
		// reset the matrix
		MultigridUtils.setValues(_matrixOfLocationsToDetach, false);
		// map the locations to detach in a matrix for faster access
		_segmentsToDetach.setLocationsToDetach(_matrixOfLocationsToDetach);
	}

	/**
	 * @param border
	 *            array of all boder points
	 * @param biomass
	 *            provides information on whether locations have biomass or not
	 * @return a set of points (segment) representing the border but having the
	 *         maximum shear value of its imediate neighbors
	 */
	private SetOfGridNodes computeTheAverageShearForEachPointInBorder(
			ArrayList border, BiomassBooleanMatrix biomass) {
		SetOfGridNodes borderWithValue = new SetOfGridNodes();
		for (Iterator iter = border.iterator(); iter.hasNext();) {
			DiscreteCoordinate p = (DiscreteCoordinate) iter.next();
			// Cycle the neighboring locations (in 2D) and compute the
			// average shear of neighbors locations in liquid
			float averageShear = 0;
			float numberOfNeighborsInLiquid = 0;
			for (int i = -1; i <= 1; i++) {
				for (int j = -1; j <= 1; j++) {
					// visit nighbors and check if they are empty and
					// add the shear value for averaging
					if (!biomass.isCarrierOrBiomass2D(p.i + i, p.j + j)) {
						averageShear += TubeVelocityField.getInstance()
								.getShearStress().getValueAt(p.i + i, p.j + j,
										1);
						numberOfNeighborsInLiquid++;
					}
				}
			}
			averageShear /= numberOfNeighborsInLiquid;
			// create the position with associated shear value
			borderWithValue.add(new GridNode2DWithValue(p, averageShear));
		}
		return borderWithValue;
	}

	/**
	 * looks in the matrix if this bacteria is located in a region to detach
	 * 
	 * @param b
	 * @return
	 */
	public boolean checkIfBacteriaShouldDetach(BiomassParticle b) {
		int i = MultigridVariable.snapToDiscreteI(b.getCenterX());
		int j = MultigridVariable.snapToDiscreteJ(b.getCenterY());
		return _matrixOfLocationsToDetach[i][j][1];
	}

	/**
	 * Write the border with value to disk
	 * 
	 * @param dir
	 */
	public void writeBorderWithValue(String dir) {
		if (_borderWithValue != null)
			writeObjectToFile(_borderWithValue, dir + "borderWithValue.txt");
	}

	/**
	 * Write each of the segments to detach (if any) to disk. each segment will
	 * be writen on a separate file.
	 * 
	 * @param dir
	 */
	public void writeSegmentsToDetach(String dir) {
		if (_segmentsToDetach != null)
			for (int i = 0; i < _segmentsToDetach.size(); i++) {
				SegmentToDetach s = (SegmentToDetach) (_segmentsToDetach.get(i));
				String fn = "segment2detach_" + (i + 1) + "_radius_" + s.radius
						+ ".txt";
				writeObjectToFile(s, dir + fn);
				String fn2 = "areas2detach_" + (i + 1) + "_radius_" + s.radius
						+ ".txt";
				writeObjectToFile(s._locationsToDetach, dir + fn2);
			}
	}

	/**
	 * Write the biomass matrix to file
	 * 
	 * @param dir
	 */
	public void writeBiomassMatrix(String dir) {
		if (_biomass != null)
			writeObjectToFile(_biomass, dir + "biomassWithCarrier.txt");
	}

	/**
	 * Auxiliary method for debugging TODO delete after testing
	 * 
	 * @param o
	 * @param fn
	 */
	private void writeObjectToFile(Object o, String fn) {
		// create the file and write matrix to it
		File f = new File(fn);
		try {
			java.io.FileWriter fr = new java.io.FileWriter(f);
			fr.write(o.toString());
			fr.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public boolean detachmentIsOff() {
		return _breakingShearStress == 0;
	}

}