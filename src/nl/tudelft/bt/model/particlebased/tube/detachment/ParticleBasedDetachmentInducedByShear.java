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
import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.multigrid.TubeVelocityField;
import nl.tudelft.bt.model.particlebased.tube.TubeBiomassParticleContainer;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * Includes the algorithms for determining which pices of biofilm will detach as
 * a consequence of the shear stress in the tube
 * 
 * @author jxavier
 */
public class ParticleBasedDetachmentInducedByShear implements DetachmentHandler {
	private ArrayList _border;

	private SetOfGridNodes _borderWithValue;

	private Segments _segmentsToDetach;

	private BiomassBooleanMatrix _biomass;

	private TubeBiomassParticleContainer _tube;

	private float _breakingShearStress;

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
			return (theta - n.theta + ExtraMath.PI2) % ExtraMath.PI2;
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

		float _shear;

		float _radius;

		float _shearPathIntegral;

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
			_radius = r;
		}

		/**
		 * Constructor that only initializes the value of the radius
		 * 
		 * @param r
		 */
		public SetOfGridNodes(float r) {
			_radius = r;
		}

		/**
		 * set the value of the radius for this segment
		 * 
		 * @param r
		 */
		void setRadius(float r) {
			_radius = r;
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
			_shear = 0;
			_radius = 0;
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
		 * Update the value of _shearPathIntegral with the integral for the
		 * whole segment
		 */
		void preDetermineLineIntegral() {
			_shearPathIntegral = computeShearPathIntegral(0, size() - 1);
		}

		/**
		 * Use the predertmined value of the line integral to compute the
		 * average shear as a function of the radius
		 * 
		 * @param r
		 * @return the average shear
		 */
		float computeShearOfCrown(float r) {
			float archDistance = ExtraMath.PI2 * r;
			// finish computing the average, dividing by the path length
			return _shearPathIntegral / archDistance;
		}

		/**
		 * Compute the average shear applied on the arch defined by the entire
		 * segment
		 */
		void computeShear() {
			int iLow = 0;
			int iHigh = size() - 1;
			GridNode2DWithValue low = get(iLow);
			GridNode2DWithValue high = get(iHigh);
			// determine the angle
			float angle = high.rotation(low);
			// compute the arch distance between the two points
			float archDistance = angle * _radius;
			// finish computing the average, dividing by the path length
			_shear = computeShearPathIntegral(iLow, iHigh) / archDistance;
		}

		/**
		 * @param iLow
		 * @param iHigh
		 * @return the path integral
		 */
		private float computeShearPathIntegral(int iLow, int iHigh) {
			float shear = 0;
			int nPoints = this.size();
			int i = iLow;
			// compute the line path
			while (i != iHigh) {
				int iNext = (i + 1) % nPoints;
				float d = get(i).distance(get(iNext));
				shear += get(i).averageShear(get(iNext)) * d;
				i = iNext;
			}
			return shear;
		}

		/**
		 * The average shear applien on an acrh beteween points iLow and iHigh
		 * 
		 * @param iLow
		 * @param iHigh
		 * @param r
		 * @return the average shear
		 */
		float computeAverageArchShear(int iLow, int iHigh, float r) {
			GridNode2DWithValue low = get(iLow);
			GridNode2DWithValue high = get(iHigh);
			// determine the angle
			float angle = low.rotation(high);
			// compute the arch distance between the two points
			float archDistance = angle * r;
			// finish computing the average, dividing by the path length
			return computeShearPathIntegral(iLow, iHigh) / archDistance;
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
			return searchFor(i, j) > -1;
		}

		/**
		 * @param i
		 * @param j
		 * @return the index of point with coodrinates i,j or -1 if the point is
		 *         not found
		 */
		public int searchFor(int i, int j) {
			for (int n = 0; n < this.size(); n++) {
				if (this.get(n).isSameLocation(i, j))
					return n;
			}
			return -1;
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
		 * Initialize a segment to detach with the date from a set that was
		 * evaluated
		 * 
		 * @param s
		 */
		public SegmentToDetach(SetOfGridNodes s) {
			super(s, s._radius);
			// find the locations to detach
			_locationsToDetach = new SetOfGridNodes(_radius);
			SetOfGridNodes trial = new SetOfGridNodes(_radius);
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
							// was
							// not found before
							if (_biomass.isCarrierOrBiomass2D(iPoint, jPoint)
									& (!_locationsToDetach.contains(iPoint,
											jPoint))
									& (!trial.contains(iPoint, jPoint))
									& (!found.contains(iPoint, jPoint))
									& isInsideRadius(_radius, iPoint, jPoint)) {
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
			if (s._shear > _breakingShearStress)
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
	}

	/**
	 * Initialize the detachment handler
	 */
	public ParticleBasedDetachmentInducedByShear(float detachmentShearStress) {
		_tube = (TubeBiomassParticleContainer) (Model.model().biomassContainer);
		_breakingShearStress = detachmentShearStress;
	}

	/**
	 * Updates the references to border points and the average shear acting on
	 * each point. Border is defined as a set of nodes of a grid (with the
	 * resolution of the finest grid in the multigrid) located in the border of
	 * the biofilm
	 */
	public void updateShearOfBorderPoints() {
		// get the updated matrix
		_biomass = Model.model().getBiomassBooleanMatrix();
		_border = _biomass.findTheBorderPointsInTube();
		// Compute the shear in each point;
		_borderWithValue = computeTheAverageShearForEachPointInBorder(_border,
				_biomass);
		_borderWithValue.preDetermineLineIntegral();
	}

	/**
	 * Compute the shear either by checking the value (if point is at border) or
	 * by computing the average shear of an arch passing though the point
	 * 
	 * @param x
	 * @param y
	 * @return the shear at that point
	 */
	public float computeShear(float x, float y) {
		int i = MultigridVariable.snapToDiscreteI(x);
		int j = MultigridVariable.snapToDiscreteI(y);
		float r = _tube.rConvertToPolar(x, y);
		// if location is in border (index > 0),
		// use the border point for the search
		int index = _borderWithValue.searchFor(i, j);
		if (index > 0) {
			float rGridNode = _borderWithValue.get(index).r;
			//
			int nPoints = _borderWithValue.size();
			// Check if point is on the top of finger like structures
			float nextRTop = _borderWithValue.get((index + 1) % nPoints).r;
			float nextRBottom = _borderWithValue.get((index - 1 + nPoints)
					% nPoints).r;
			if ((nextRTop >= rGridNode) & (nextRBottom >= rGridNode)) {
				return _borderWithValue.get(index).value;
			}
			// otherwise proceed
			return averageArchShear2(index, rGridNode);
		}
		// if it is not compute the average shear applied on the arch
		// Find the position in the border with the closest angle to this
		// position that is inside the radius
		float minDistance = Float.POSITIVE_INFINITY;
		int iClosest = -1;
		boolean imersed = true;
		for (int n = 0; n < _borderWithValue.size(); n++) {
			GridNode2DWithValue g = _borderWithValue.get(n);
			float distance = ExtraMath.pointDistance(x, y, g.x, g.y);
			if (g.r < r) {
				if (distance < minDistance) {
					minDistance = distance;
					iClosest = n;
				}
			} else
				imersed = false;
		}
		if ((imersed) | (iClosest == -1))
			return _borderWithValue.computeShearOfCrown(r);
		return averageArchShear2(iClosest, r);
	}

	/**
	 * @param iClosest
	 * @param r
	 * @return
	 */
	public float averageArchShear2(int iClosest, float r) {
		int nPoints = _borderWithValue.size();
		int iHigh = iClosest;
		while (_borderWithValue.get(iHigh).r <= r) {
			iHigh = (iHigh + 1) % nPoints;
			if (iHigh == iClosest)
				break;
		}
		// find the lower value
		int iLow = iClosest;
		while (_borderWithValue.get(iLow).r <= r) {
			iLow = (iLow - 1 + nPoints) % nPoints;
			if (iLow == iClosest)
				break;
		}
		// case is that the particle is almost completely imersed in the
		// biofilm
		if ((iLow == iHigh)) {
			return _borderWithValue.computeShearOfCrown(r);
		}
		return _borderWithValue.computeAverageArchShear(iLow, iHigh, r);
	}

	public boolean detachmentIsOff() {
		// Detachment is not computed if value is 0 for the
		// _breakingShearStress
		return _breakingShearStress == 0;
	}

	public float averageArchShear(float x, float y) {
		float r = _tube.rConvertToPolar(x, y);
		// Find the position in the border with the closest angle to this
		// position that is inside the radius
		float minDistance = Float.POSITIVE_INFINITY;
		int iClosest = -1;
		for (int i = 0; i < _borderWithValue.size(); i++) {
			GridNode2DWithValue n = _borderWithValue.get(i);
			float distance = ExtraMath.pointDistance(x, y, n.x, n.y);
			if (n.r < r)
				if (distance < minDistance) {
					minDistance = distance;
					iClosest = i;
				}
		}
		// if no point was found, then the case is that the point is in the
		// border (on the top of finger like structures)
		if (iClosest == -1) {
			int i = MultigridVariable.snapToDiscreteI(x);
			int j = MultigridVariable.snapToDiscreteI(y);
			int index = _borderWithValue.searchFor(i, j);
			if (index > 0)
				return _borderWithValue.get(index).value;
		}
		// otherwise, proceed
		int nPoints = _borderWithValue.size();
		// find the upper value
		int iHigh = iClosest;
		while (_borderWithValue.get(iHigh).r < r) {
			iHigh = (iHigh + 1) % nPoints;
			if (iHigh == iClosest)
				break;
		}
		// find the lower value
		int iLow = iClosest;
		while (_borderWithValue.get(iLow).r < r) {
			iLow = (iLow - 1 + nPoints) % nPoints;
			if (iLow == iClosest)
				break;
		}

		// case is that the particle is completely imersed in the
		// biofilm
		if ((iLow == iHigh)) {
			return _borderWithValue.computeShearOfCrown(r);
		}
		return _borderWithValue.computeAverageArchShear(iLow, iHigh, r);
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
		for (int i = 0; i < _segmentsToDetach.size(); i++) {
			SegmentToDetach s = (SegmentToDetach) (_segmentsToDetach.get(i));
			String fn = "segment2detach_" + (i + 1) + "_radius_" + s._radius
					+ ".txt";
			writeObjectToFile(s, dir + fn);
			String fn2 = "areas2detach_" + (i + 1) + "_radius_" + s._radius
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

	/**
	 * @return the value of the breaking shear stress
	 */
	public float getBreakingShearStress() {
		return _breakingShearStress;
	}

}