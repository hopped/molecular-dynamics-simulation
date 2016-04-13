/*
 * Copyright (c) 2014      Christoph Niethammer <christoph.niethammer@gmail.com>
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER
 */

#ifndef OBJECTS_H
#define OBJECTS_H

#ifdef __cplusplus
class Object {
public:
	virtual ~Object(){}
	/** Determines if the given point is inside the object */
	virtual bool isInside(double r[3]) = 0;
	/** Get lower corner of a bounding box around the object */
	virtual void getBboxMin(double rmin[3]) = 0;
	/** Get upper corner of a bounding box around the object */
	virtual void getBboxMax(double rmax[3]) = 0;
};

/** Class implementing a cuboid */
class Cuboid : public Object {
public:
	Cuboid(double lower[3], double upper[3]) {
		for(int d = 0; d < 3; d++) {
			_lowerCorner[d] = lower[d];
			_upperCorner[d] = upper[d];
		}
	}

	bool isInside(double r[3]) {
		return (_lowerCorner[0] <= r[0] && r[0] < _upperCorner[0])
			&& (_lowerCorner[1] <= r[1] && r[1] < _upperCorner[1])
			&& (_lowerCorner[2] <= r[2] && r[2] < _upperCorner[2]);
	}

	void getBboxMin(double rmin[3]) {
		for(int d = 0; d < 3; d++) {
			rmin[d] = _lowerCorner[d];
		}
	}
	void getBboxMax(double rmax[3]) {
		for(int d = 0; d < 3; d++) {
			rmax[d] = _upperCorner[d];
		}
	}

private:
	double _lowerCorner[3];
	double _upperCorner[3];
};

/** Class implementing a sphere */
class Sphere : public Object {
public:
	Sphere(double center[3], double r) : _radius(r), _radiusSquare(r*r) {
		for(int d = 0; d < 3; d++) {
			_center[d] = center[d];
		}
	}

	bool isInside(double r[3]) {
		double dr2 = 0.0;
		double dr[3];
		for(int d = 0; d < 3; d++) {
			dr[d] = r[d] - _center[d];
			dr2 += dr[d] * dr[d];
		}
		return (dr2 < _radiusSquare);
	}

	void getBboxMin(double rmin[3]) {
		for(int d = 0; d < 3; d++) {
			rmin[d] = _center[d] - _radius;
		}
	}
	void getBboxMax(double rmax[3]) {
		for(int d = 0; d < 3; d++) {
			rmax[d] = _center[d] + _radius;
		}
	}

private:
	double _center[3];
	double _radius;
	double _radiusSquare;
};

/** Class implementing a cyliner */
class Cylinder : public Object {
public:
	/** Constructor
	 * @param[in]  centerBase   Center of the circle of the lower base of the cylinder.
	 * @param[in]  radius       Raius of the cylinder (x-y-axis)
	 * @param[in]  height       Height of the cylinder (z-axis)
	 */
	Cylinder(double centerBase[3], double radius, double height) : _radius(radius), _height(height) {
		for(int d = 0; d < 3; d++) {
			_centerBase[d] = centerBase[d];
		}
	}

	bool isInside(double r[3]) {
		double dr[2];
		dr[0] = r[0] - _centerBase[0];
		dr[1] = r[1] - _centerBase[1];
		return (dr[0]*dr[0] + dr[1]*dr[1] < _radius * _radius) && (r[2] >= _centerBase[0]) && (r[2] < _centerBase[2] + _height);
	}

	void getBboxMin(double rmin[3]) {
		rmin[0] = _centerBase[0] - _radius;
		rmin[1] = _centerBase[1] - _radius;
		rmin[2] = _centerBase[2];
	}
	void getBboxMax(double rmax[3]) {
		rmax[0] = _centerBase[0] + _radius;
		rmax[1] = _centerBase[1] + _radius;
		rmax[2] = _centerBase[2] + _height;
	}

private:
	double _radius;
	double _height;
	double _centerBase[3];
};

/** Abstract class to combine two objects */
class ObjectCombiner : public Object {
public:
	ObjectCombiner(Object *ob1, Object *ob2) : _ob1(ob1), _ob2(ob2) {}

	bool isInside(double r[3]) {
		return _ob1->isInside(r) && _ob2->isInside(r);
	}

	void getBboxMin(double rmin[3]) {
		double rmin1[3], rmin2[3];
		_ob1->getBboxMin(rmin1);
		_ob2->getBboxMin(rmin2);
		for(int d = 0; d < 3; d++) {
			rmin[d] = (rmin1[d] < rmin2[d]) ? rmin1[d] : rmin2[d] ;
		}
	}
	void getBboxMax(double rmax[3]) {
		double rmax1[3], rmax2[3];
		_ob1->getBboxMax(rmax1);
		_ob2->getBboxMax(rmax2);
		for(int d = 0; d < 3; d++) {
			rmax[d] = (rmax1[d] > rmax2[d]) ? rmax1[d] : rmax2[d] ;
		}
	}

private:
	Object* _ob1;
	Object* _ob2;
};

/** Abstract class to subtract one object from another */
class ObjectSubtractor : public Object {
public:
	/** Constructor
	 * @param[in]  original_ob  The original object.
	 * @param[in]  subtract_ob  The object which shall be subtract from the original object.
	 */
	ObjectSubtractor(Object *original_ob, Object *subtract_ob) : _ob1(original_ob), _ob2(subtract_ob) {}

	bool isInside(double r[3]) {
		return _ob1->isInside(r) && (!_ob2->isInside(r));
	}

	void getBboxMin(double rmin[3]) {
		_ob1->getBboxMin(rmin);
	}
	void getBboxMax(double rmax[3]) {
		_ob1->getBboxMax(rmax);
	}

private:
	Object* _ob1;
	Object* _ob2;
};

#endif /* __cplusplus */

#endif /* OBJECTS_H */
