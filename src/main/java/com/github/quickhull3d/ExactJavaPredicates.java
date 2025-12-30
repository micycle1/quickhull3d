package com.github.quickhull3d;

class ExactJavaPredicates {

	// from ProGAL

	static {
		exactinit();
	}

	private ExactJavaPredicates() {
	}

	static double orient(Vertex p0, Vertex p1, Vertex p2, Vertex q) {
		return orient3d(p0.x(), p0.y(), p0.z(), p1.x(), p1.y(), p1.z(), p2.x(), p2.y(), p2.z(), q.x(), q.y(), q.z());
	}

	static double orient3d(double pax, double pay, double paz, double pbx, double pby, double pbz, double pcx, double pcy, double pcz, double pdx, double pdy,
			double pdz) {
		double adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz;
		double bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady;
		double det;
		double permanent, errbound;

		adx = pax - pdx;
		bdx = pbx - pdx;
		cdx = pcx - pdx;
		ady = pay - pdy;
		bdy = pby - pdy;
		cdy = pcy - pdy;
		adz = paz - pdz;
		bdz = pbz - pdz;
		cdz = pcz - pdz;

		bdxcdy = bdx * cdy;
		cdxbdy = cdx * bdy;

		cdxady = cdx * ady;
		adxcdy = adx * cdy;

		adxbdy = adx * bdy;
		bdxady = bdx * ady;

		det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);

		permanent = (((bdxcdy) >= 0.0 ? (bdxcdy) : -(bdxcdy)) + ((cdxbdy) >= 0.0 ? (cdxbdy) : -(cdxbdy))) * ((adz) >= 0.0 ? (adz) : -(adz))
				+ (((cdxady) >= 0.0 ? (cdxady) : -(cdxady)) + ((adxcdy) >= 0.0 ? (adxcdy) : -(adxcdy))) * ((bdz) >= 0.0 ? (bdz) : -(bdz))
				+ (((adxbdy) >= 0.0 ? (adxbdy) : -(adxbdy)) + ((bdxady) >= 0.0 ? (bdxady) : -(bdxady))) * ((cdz) >= 0.0 ? (cdz) : -(cdz));
		errbound = o3derrboundA * permanent;

		if ((det > errbound) || (-det > errbound)) {
//			System.out.println("Returning early " + det + " .. " + errbound);
			return det;
		}
		// Only allocate arrays when adaptive precision is needed (rare path)
		return orient3dadaptHelper(pax, pay, paz, pbx, pby, pbz, pcx, pcy, pcz, pdx, pdy, pdz, permanent);
	}

	private static double orient3dadaptHelper(double pax, double pay, double paz, double pbx, double pby, double pbz, double pcx, double pcy, double pcz,
			double pdx, double pdy, double pdz, double permanent) {
		double[] pa = new double[3];
		double[] pb = new double[3];
		double[] pc = new double[3];
		double[] pd = new double[3];

		pa[0] = pax;
		pa[1] = pay;
		pa[2] = paz;
		pb[0] = pbx;
		pb[1] = pby;
		pb[2] = pbz;
		pc[0] = pcx;
		pc[1] = pcy;
		pc[2] = pcz;
		pd[0] = pdx;
		pd[1] = pdy;
		pd[2] = pdz;

		return orient3dadapt(pa, pb, pc, pd, permanent);
	}

	// Fast version without error bounds checking
	static double orient3dfast(double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz, double dx, double dy,
			double dz) {
		double adx = ax - dx;
		double bdx = bx - dx;
		double cdx = cx - dx;
		double ady = ay - dy;
		double bdy = by - dy;
		double cdy = cy - dy;
		double adz = az - dz;
		double bdz = bz - dz;
		double cdz = cz - dz;

		return adx * (bdy * cdz - bdz * cdy) + bdx * (cdy * adz - cdz * ady) + cdx * (ady * bdz - adz * bdy);
	}

	private static double splitter;
	private static double epsilon;

	private static double resulterrbound;
	private static double o3derrboundA, o3derrboundB, o3derrboundC;

	private static void exactinit() {
		double half;
		double check, lastcheck;
		int every_other;

		every_other = 1;
		half = 0.5;
		epsilon = 1.0;
		splitter = 1.0;
		check = 1.0;

		do {
			lastcheck = check;
			epsilon *= half;
			if (every_other != 0) {
				splitter *= 2.0;
			}
			every_other = every_other == 0 ? 1 : 0;
			check = 1.0 + epsilon;
		} while ((check != 1.0) && (check != lastcheck));
		splitter += 1.0;

		resulterrbound = (3.0 + 8.0 * epsilon) * epsilon;
		o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon;
		o3derrboundB = (3.0 + 28.0 * epsilon) * epsilon;
		o3derrboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon;
	}

	private static int fast_expansion_sum_zeroelim(int elen, double[] e, int flen, double[] f, double[] h) {
		double Q;
		double Qnew;
		double hh;
		double bvirt;
		double avirt, bround, around;
		int eindex, findex, hindex;
		double enow, fnow;

		enow = e[0];
		fnow = f[0];
		eindex = findex = 0;
		if ((fnow > enow) == (fnow > -enow)) {
			Q = enow;
			enow = e[++eindex];
		} else {
			Q = fnow;
			fnow = f[++findex];
		}
		hindex = 0;
		if ((eindex < elen) && (findex < flen)) {
			if ((fnow > enow) == (fnow > -enow)) {
				Qnew = (double) (enow + Q);
				bvirt = Qnew - enow;
				hh = Q - bvirt;
				enow = e[++eindex];
			} else {
				Qnew = (double) (fnow + Q);
				bvirt = Qnew - fnow;
				hh = Q - bvirt;
				fnow = f[++findex];
			}
			Q = Qnew;
			if (hh != 0.0) {
				h[hindex++] = hh;
			}
			while ((eindex < elen) && (findex < flen)) {
				if ((fnow > enow) == (fnow > -enow)) {
					Qnew = (double) (Q + enow);
					bvirt = (double) (Qnew - Q);
					avirt = Qnew - bvirt;
					bround = enow - bvirt;
					around = Q - avirt;
					hh = around + bround;
					eindex++;
					if (eindex < elen) { // Check before reading
						enow = e[eindex];
					}
				} else {
					Qnew = (double) (Q + fnow);
					bvirt = (double) (Qnew - Q);
					avirt = Qnew - bvirt;
					bround = fnow - bvirt;
					around = Q - avirt;
					hh = around + bround;
					findex++;
					if (findex < flen) { // Check before reading
						fnow = f[findex];
					}
				}
				Q = Qnew;
				if (hh != 0.0) {
					h[hindex++] = hh;
				}
			}
		}
		while (eindex < elen) {
			Qnew = (double) (Q + enow);
			bvirt = (double) (Qnew - Q);
			avirt = Qnew - bvirt;
			bround = enow - bvirt;
			around = Q - avirt;
			hh = around + bround;
			Q = Qnew;
			if (hh != 0.0) {
				h[hindex++] = hh;
			}
			eindex++;
			if (eindex < elen) { // Check before reading
				enow = e[eindex];
			}
		}

		while (findex < flen) {
			Qnew = (double) (Q + fnow);
			bvirt = (double) (Qnew - Q);
			avirt = Qnew - bvirt;
			bround = fnow - bvirt;
			around = Q - avirt;
			hh = around + bround;
			Q = Qnew;
			if (hh != 0.0) {
				h[hindex++] = hh;
			}
			findex++;
			if (findex < flen) { // Check before reading
				fnow = f[findex];
			}
		}
		if ((Q != 0.0) || (hindex == 0)) {
			h[hindex++] = Q;
		}
		return hindex;
	}

	private static int scale_expansion_zeroelim(int elen, double[] e, double b, double[] h) {
		double Q, sum;
		double hh;
		double product1;
		double product0;
		int eindex, hindex;
		double enow;
		double bvirt;
		double avirt, bround, around;
		double c;
		double abig;
		double ahi, alo, bhi, blo;
		double err1, err2, err3;

		c = (double) (splitter * b);
		abig = (double) (c - b);
		bhi = c - abig;
		blo = b - bhi;
		Q = (double) (e[0] * b);
		c = (double) (splitter * e[0]);
		abig = (double) (c - e[0]);
		ahi = c - abig;
		alo = e[0] - ahi;
		err1 = Q - (ahi * bhi);
		err2 = err1 - (alo * bhi);
		err3 = err2 - (ahi * blo);
		hh = (alo * blo) - err3;
		hindex = 0;
		if (hh != 0) {
			h[hindex++] = hh;
		}
		for (eindex = 1; eindex < elen; eindex++) {
			enow = e[eindex];
			product1 = (double) (enow * b);
			c = (double) (splitter * enow);
			abig = (double) (c - enow);
			ahi = c - abig;
			alo = enow - ahi;
			err1 = product1 - (ahi * bhi);
			err2 = err1 - (alo * bhi);
			err3 = err2 - (ahi * blo);
			product0 = (alo * blo) - err3;
			sum = (double) (Q + product0);
			bvirt = (double) (sum - Q);
			avirt = sum - bvirt;
			bround = product0 - bvirt;
			around = Q - avirt;
			hh = around + bround;
			if (hh != 0) {
				h[hindex++] = hh;
			}
			Q = (double) (product1 + sum);
			bvirt = Q - product1;
			hh = sum - bvirt;
			if (hh != 0) {
				h[hindex++] = hh;
			}
		}
		if ((Q != 0.0) || (hindex == 0)) {
			h[hindex++] = Q;
		}
		return hindex;
	}

	private static double estimate(int elen, double[] e) {
		double Q;
		int eindex;

		Q = e[0];
		for (eindex = 1; eindex < elen; eindex++) {
			Q += e[eindex];
		}
		return Q;
	}

	private static double orient3dadapt(double[] pa, double[] pb, double[] pc, double[] pd, double permanent) {
		double adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz;
		double det, errbound;

		double bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
		double bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
		double[] bc = new double[4], ca = new double[4], ab = new double[4];
		double bc3, ca3, ab3;
		double[] adet = new double[8], bdet = new double[8], cdet = new double[8];
		int alen, blen, clen;
		double[] abdet = new double[16];
		int ablen;
		double[] finnow, finother, finswap;
		double[] fin1 = new double[192], fin2 = new double[192];
		int finlength;

		double adxtail, bdxtail, cdxtail;
		double adytail, bdytail, cdytail;
		double adztail, bdztail, cdztail;
		double at_blarge, at_clarge;
		double bt_clarge, bt_alarge;
		double ct_alarge, ct_blarge;
		double[] at_b = new double[4], at_c = new double[4], bt_c = new double[4], bt_a = new double[4], ct_a = new double[4], ct_b = new double[4];
		int at_blen, at_clen, bt_clen, bt_alen, ct_alen, ct_blen;
		double bdxt_cdy1, cdxt_bdy1, cdxt_ady1;
		double adxt_cdy1, adxt_bdy1, bdxt_ady1;
		double bdxt_cdy0, cdxt_bdy0, cdxt_ady0;
		double adxt_cdy0, adxt_bdy0, bdxt_ady0;
		double bdyt_cdx1, cdyt_bdx1, cdyt_adx1;
		double adyt_cdx1, adyt_bdx1, bdyt_adx1;
		double bdyt_cdx0, cdyt_bdx0, cdyt_adx0;
		double adyt_cdx0, adyt_bdx0, bdyt_adx0;
		double[] bct = new double[8], cat = new double[8], abt = new double[8];
		int bctlen, catlen, abtlen;
		double bdxt_cdyt1, cdxt_bdyt1, cdxt_adyt1;
		double adxt_cdyt1, adxt_bdyt1, bdxt_adyt1;
		double bdxt_cdyt0, cdxt_bdyt0, cdxt_adyt0;
		double adxt_cdyt0, adxt_bdyt0, bdxt_adyt0;
		double[] u = new double[4], v = new double[12], w = new double[16];
		double u3;
		int vlength, wlength;
		double negate;

		double bvirt;
		double avirt, bround, around;
		double c;
		double abig;
		double ahi, alo, bhi, blo;
		double err1, err2, err3;
		double _i, _j, _k;
		double _0;

		adx = (double) (pa[0] - pd[0]);
		bdx = (double) (pb[0] - pd[0]);
		cdx = (double) (pc[0] - pd[0]);
		ady = (double) (pa[1] - pd[1]);
		bdy = (double) (pb[1] - pd[1]);
		cdy = (double) (pc[1] - pd[1]);
		adz = (double) (pa[2] - pd[2]);
		bdz = (double) (pb[2] - pd[2]);
		cdz = (double) (pc[2] - pd[2]);

		bdxcdy1 = (double) (bdx * cdy);
		c = (double) (splitter * bdx);
		abig = (double) (c - bdx);
		ahi = c - abig;
		alo = bdx - ahi;
		c = (double) (splitter * cdy);
		abig = (double) (c - cdy);
		bhi = c - abig;
		blo = cdy - bhi;
		err1 = bdxcdy1 - (ahi * bhi);
		err2 = err1 - (alo * bhi);
		err3 = err2 - (ahi * blo);
		bdxcdy0 = (alo * blo) - err3;
		cdxbdy1 = (double) (cdx * bdy);
		c = (double) (splitter * cdx);
		abig = (double) (c - cdx);
		ahi = c - abig;
		alo = cdx - ahi;
		c = (double) (splitter * bdy);
		abig = (double) (c - bdy);
		bhi = c - abig;
		blo = bdy - bhi;
		err1 = cdxbdy1 - (ahi * bhi);
		err2 = err1 - (alo * bhi);
		err3 = err2 - (ahi * blo);
		cdxbdy0 = (alo * blo) - err3;
		_i = (double) (bdxcdy0 - cdxbdy0);
		bvirt = (double) (bdxcdy0 - _i);
		avirt = _i + bvirt;
		bround = bvirt - cdxbdy0;
		around = bdxcdy0 - avirt;
		bc[0] = around + bround;
		_j = (double) (bdxcdy1 + _i);
		bvirt = (double) (_j - bdxcdy1);
		avirt = _j - bvirt;
		bround = _i - bvirt;
		around = bdxcdy1 - avirt;
		_0 = around + bround;
		_i = (double) (_0 - cdxbdy1);
		bvirt = (double) (_0 - _i);
		avirt = _i + bvirt;
		bround = bvirt - cdxbdy1;
		around = _0 - avirt;
		bc[1] = around + bround;
		bc3 = (double) (_j + _i);
		bvirt = (double) (bc3 - _j);
		avirt = bc3 - bvirt;
		bround = _i - bvirt;
		around = _j - avirt;
		bc[2] = around + bround;
		bc[3] = bc3;
		alen = scale_expansion_zeroelim(4, bc, adz, adet);

		cdxady1 = (double) (cdx * ady);
		c = (double) (splitter * cdx);
		abig = (double) (c - cdx);
		ahi = c - abig;
		alo = cdx - ahi;
		c = (double) (splitter * ady);
		abig = (double) (c - ady);
		bhi = c - abig;
		blo = ady - bhi;
		err1 = cdxady1 - (ahi * bhi);
		err2 = err1 - (alo * bhi);
		err3 = err2 - (ahi * blo);
		cdxady0 = (alo * blo) - err3;
		adxcdy1 = (double) (adx * cdy);
		c = (double) (splitter * adx);
		abig = (double) (c - adx);
		ahi = c - abig;
		alo = adx - ahi;
		c = (double) (splitter * cdy);
		abig = (double) (c - cdy);
		bhi = c - abig;
		blo = cdy - bhi;
		err1 = adxcdy1 - (ahi * bhi);
		err2 = err1 - (alo * bhi);
		err3 = err2 - (ahi * blo);
		adxcdy0 = (alo * blo) - err3;
		_i = (double) (cdxady0 - adxcdy0);
		bvirt = (double) (cdxady0 - _i);
		avirt = _i + bvirt;
		bround = bvirt - adxcdy0;
		around = cdxady0 - avirt;
		ca[0] = around + bround;
		_j = (double) (cdxady1 + _i);
		bvirt = (double) (_j - cdxady1);
		avirt = _j - bvirt;
		bround = _i - bvirt;
		around = cdxady1 - avirt;
		_0 = around + bround;
		_i = (double) (_0 - adxcdy1);
		bvirt = (double) (_0 - _i);
		avirt = _i + bvirt;
		bround = bvirt - adxcdy1;
		around = _0 - avirt;
		ca[1] = around + bround;
		ca3 = (double) (_j + _i);
		bvirt = (double) (ca3 - _j);
		avirt = ca3 - bvirt;
		bround = _i - bvirt;
		around = _j - avirt;
		ca[2] = around + bround;
		ca[3] = ca3;
		blen = scale_expansion_zeroelim(4, ca, bdz, bdet);

		adxbdy1 = (double) (adx * bdy);
		c = (double) (splitter * adx);
		abig = (double) (c - adx);
		ahi = c - abig;
		alo = adx - ahi;
		c = (double) (splitter * bdy);
		abig = (double) (c - bdy);
		bhi = c - abig;
		blo = bdy - bhi;
		err1 = adxbdy1 - (ahi * bhi);
		err2 = err1 - (alo * bhi);
		err3 = err2 - (ahi * blo);
		adxbdy0 = (alo * blo) - err3;
		bdxady1 = (double) (bdx * ady);
		c = (double) (splitter * bdx);
		abig = (double) (c - bdx);
		ahi = c - abig;
		alo = bdx - ahi;
		c = (double) (splitter * ady);
		abig = (double) (c - ady);
		bhi = c - abig;
		blo = ady - bhi;
		err1 = bdxady1 - (ahi * bhi);
		err2 = err1 - (alo * bhi);
		err3 = err2 - (ahi * blo);
		bdxady0 = (alo * blo) - err3;
		_i = (double) (adxbdy0 - bdxady0);
		bvirt = (double) (adxbdy0 - _i);
		avirt = _i + bvirt;
		bround = bvirt - bdxady0;
		around = adxbdy0 - avirt;
		ab[0] = around + bround;
		_j = (double) (adxbdy1 + _i);
		bvirt = (double) (_j - adxbdy1);
		avirt = _j - bvirt;
		bround = _i - bvirt;
		around = adxbdy1 - avirt;
		_0 = around + bround;
		_i = (double) (_0 - bdxady1);
		bvirt = (double) (_0 - _i);
		avirt = _i + bvirt;
		bround = bvirt - bdxady1;
		around = _0 - avirt;
		ab[1] = around + bround;
		ab3 = (double) (_j + _i);
		bvirt = (double) (ab3 - _j);
		avirt = ab3 - bvirt;
		bround = _i - bvirt;
		around = _j - avirt;
		ab[2] = around + bround;
		ab[3] = ab3;
		clen = scale_expansion_zeroelim(4, ab, cdz, cdet);

		ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
		finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1);

		det = estimate(finlength, fin1);
		errbound = o3derrboundB * permanent;
		if ((det >= errbound) || (-det >= errbound)) {
			return det;
		}

		bvirt = (double) (pa[0] - adx);
		avirt = adx + bvirt;
		bround = bvirt - pd[0];
		around = pa[0] - avirt;
		adxtail = around + bround;
		bvirt = (double) (pb[0] - bdx);
		avirt = bdx + bvirt;
		bround = bvirt - pd[0];
		around = pb[0] - avirt;
		bdxtail = around + bround;
		bvirt = (double) (pc[0] - cdx);
		avirt = cdx + bvirt;
		bround = bvirt - pd[0];
		around = pc[0] - avirt;
		cdxtail = around + bround;
		bvirt = (double) (pa[1] - ady);
		avirt = ady + bvirt;
		bround = bvirt - pd[1];
		around = pa[1] - avirt;
		adytail = around + bround;
		bvirt = (double) (pb[1] - bdy);
		avirt = bdy + bvirt;
		bround = bvirt - pd[1];
		around = pb[1] - avirt;
		bdytail = around + bround;
		bvirt = (double) (pc[1] - cdy);
		avirt = cdy + bvirt;
		bround = bvirt - pd[1];
		around = pc[1] - avirt;
		cdytail = around + bround;
		bvirt = (double) (pa[2] - adz);
		avirt = adz + bvirt;
		bround = bvirt - pd[2];
		around = pa[2] - avirt;
		adztail = around + bround;
		bvirt = (double) (pb[2] - bdz);
		avirt = bdz + bvirt;
		bround = bvirt - pd[2];
		around = pb[2] - avirt;
		bdztail = around + bround;
		bvirt = (double) (pc[2] - cdz);
		avirt = cdz + bvirt;
		bround = bvirt - pd[2];
		around = pc[2] - avirt;
		cdztail = around + bround;

		if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0) && (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0) && (adztail == 0.0)
				&& (bdztail == 0.0) && (cdztail == 0.0)) {
			return det;
		}

		errbound = o3derrboundC * permanent + resulterrbound * ((det) >= 0.0 ? (det) : -(det));
		det += (adz * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail)) + adztail * (bdx * cdy - bdy * cdx))
				+ (bdz * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail)) + bdztail * (cdx * ady - cdy * adx))
				+ (cdz * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail)) + cdztail * (adx * bdy - ady * bdx));
		if ((det >= errbound) || (-det >= errbound)) {
			return det;
		}

		finnow = fin1;
		finother = fin2;

		if (adxtail == 0.0) {
			if (adytail == 0.0) {
				at_b[0] = 0.0;
				at_blen = 1;
				at_c[0] = 0.0;
				at_clen = 1;
			} else {
				negate = -adytail;
				at_blarge = (double) (negate * bdx);
				c = (double) (splitter * negate);
				abig = (double) (c - negate);
				ahi = c - abig;
				alo = negate - ahi;
				c = (double) (splitter * bdx);
				abig = (double) (c - bdx);
				bhi = c - abig;
				blo = bdx - bhi;
				err1 = at_blarge - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				at_b[0] = (alo * blo) - err3;
				at_b[1] = at_blarge;
				at_blen = 2;
				at_clarge = (double) (adytail * cdx);
				c = (double) (splitter * adytail);
				abig = (double) (c - adytail);
				ahi = c - abig;
				alo = adytail - ahi;
				c = (double) (splitter * cdx);
				abig = (double) (c - cdx);
				bhi = c - abig;
				blo = cdx - bhi;
				err1 = at_clarge - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				at_c[0] = (alo * blo) - err3;
				at_c[1] = at_clarge;
				at_clen = 2;
			}
		} else {
			if (adytail == 0.0) {
				at_blarge = (double) (adxtail * bdy);
				c = (double) (splitter * adxtail);
				abig = (double) (c - adxtail);
				ahi = c - abig;
				alo = adxtail - ahi;
				c = (double) (splitter * bdy);
				abig = (double) (c - bdy);
				bhi = c - abig;
				blo = bdy - bhi;
				err1 = at_blarge - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				at_b[0] = (alo * blo) - err3;
				at_b[1] = at_blarge;
				at_blen = 2;
				negate = -adxtail;
				at_clarge = (double) (negate * cdy);
				c = (double) (splitter * negate);
				abig = (double) (c - negate);
				ahi = c - abig;
				alo = negate - ahi;
				c = (double) (splitter * cdy);
				abig = (double) (c - cdy);
				bhi = c - abig;
				blo = cdy - bhi;
				err1 = at_clarge - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				at_c[0] = (alo * blo) - err3;
				at_c[1] = at_clarge;
				at_clen = 2;
			} else {
				adxt_bdy1 = (double) (adxtail * bdy);
				c = (double) (splitter * adxtail);
				abig = (double) (c - adxtail);
				ahi = c - abig;
				alo = adxtail - ahi;
				c = (double) (splitter * bdy);
				abig = (double) (c - bdy);
				bhi = c - abig;
				blo = bdy - bhi;
				err1 = adxt_bdy1 - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				adxt_bdy0 = (alo * blo) - err3;
				adyt_bdx1 = (double) (adytail * bdx);
				c = (double) (splitter * adytail);
				abig = (double) (c - adytail);
				ahi = c - abig;
				alo = adytail - ahi;
				c = (double) (splitter * bdx);
				abig = (double) (c - bdx);
				bhi = c - abig;
				blo = bdx - bhi;
				err1 = adyt_bdx1 - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				adyt_bdx0 = (alo * blo) - err3;
				_i = (double) (adxt_bdy0 - adyt_bdx0);
				bvirt = (double) (adxt_bdy0 - _i);
				avirt = _i + bvirt;
				bround = bvirt - adyt_bdx0;
				around = adxt_bdy0 - avirt;
				at_b[0] = around + bround;
				_j = (double) (adxt_bdy1 + _i);
				bvirt = (double) (_j - adxt_bdy1);
				avirt = _j - bvirt;
				bround = _i - bvirt;
				around = adxt_bdy1 - avirt;
				_0 = around + bround;
				_i = (double) (_0 - adyt_bdx1);
				bvirt = (double) (_0 - _i);
				avirt = _i + bvirt;
				bround = bvirt - adyt_bdx1;
				around = _0 - avirt;
				at_b[1] = around + bround;
				at_blarge = (double) (_j + _i);
				bvirt = (double) (at_blarge - _j);
				avirt = at_blarge - bvirt;
				bround = _i - bvirt;
				around = _j - avirt;
				at_b[2] = around + bround;

				at_b[3] = at_blarge;
				at_blen = 4;
				adyt_cdx1 = (double) (adytail * cdx);
				c = (double) (splitter * adytail);
				abig = (double) (c - adytail);
				ahi = c - abig;
				alo = adytail - ahi;
				c = (double) (splitter * cdx);
				abig = (double) (c - cdx);
				bhi = c - abig;
				blo = cdx - bhi;
				err1 = adyt_cdx1 - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				adyt_cdx0 = (alo * blo) - err3;
				adxt_cdy1 = (double) (adxtail * cdy);
				c = (double) (splitter * adxtail);
				abig = (double) (c - adxtail);
				ahi = c - abig;
				alo = adxtail - ahi;
				c = (double) (splitter * cdy);
				abig = (double) (c - cdy);
				bhi = c - abig;
				blo = cdy - bhi;
				err1 = adxt_cdy1 - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				adxt_cdy0 = (alo * blo) - err3;
				_i = (double) (adyt_cdx0 - adxt_cdy0);
				bvirt = (double) (adyt_cdx0 - _i);
				avirt = _i + bvirt;
				bround = bvirt - adxt_cdy0;
				around = adyt_cdx0 - avirt;
				at_c[0] = around + bround;
				_j = (double) (adyt_cdx1 + _i);
				bvirt = (double) (_j - adyt_cdx1);
				avirt = _j - bvirt;
				bround = _i - bvirt;
				around = adyt_cdx1 - avirt;
				_0 = around + bround;
				_i = (double) (_0 - adxt_cdy1);
				bvirt = (double) (_0 - _i);
				avirt = _i + bvirt;
				bround = bvirt - adxt_cdy1;
				around = _0 - avirt;
				at_c[1] = around + bround;
				at_clarge = (double) (_j + _i);
				bvirt = (double) (at_clarge - _j);
				avirt = at_clarge - bvirt;
				bround = _i - bvirt;
				around = _j - avirt;
				at_c[2] = around + bround;

				at_c[3] = at_clarge;
				at_clen = 4;
			}
		}
		if (bdxtail == 0.0) {
			if (bdytail == 0.0) {
				bt_c[0] = 0.0;
				bt_clen = 1;
				bt_a[0] = 0.0;
				bt_alen = 1;
			} else {
				negate = -bdytail;
				bt_clarge = (double) (negate * cdx);
				c = (double) (splitter * negate);
				abig = (double) (c - negate);
				ahi = c - abig;
				alo = negate - ahi;
				c = (double) (splitter * cdx);
				abig = (double) (c - cdx);
				bhi = c - abig;
				blo = cdx - bhi;
				err1 = bt_clarge - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				bt_c[0] = (alo * blo) - err3;
				bt_c[1] = bt_clarge;
				bt_clen = 2;
				bt_alarge = (double) (bdytail * adx);
				c = (double) (splitter * bdytail);
				abig = (double) (c - bdytail);
				ahi = c - abig;
				alo = bdytail - ahi;
				c = (double) (splitter * adx);
				abig = (double) (c - adx);
				bhi = c - abig;
				blo = adx - bhi;
				err1 = bt_alarge - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				bt_a[0] = (alo * blo) - err3;
				bt_a[1] = bt_alarge;
				bt_alen = 2;
			}
		} else {
			if (bdytail == 0.0) {
				bt_clarge = (double) (bdxtail * cdy);
				c = (double) (splitter * bdxtail);
				abig = (double) (c - bdxtail);
				ahi = c - abig;
				alo = bdxtail - ahi;
				c = (double) (splitter * cdy);
				abig = (double) (c - cdy);
				bhi = c - abig;
				blo = cdy - bhi;
				err1 = bt_clarge - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				bt_c[0] = (alo * blo) - err3;
				bt_c[1] = bt_clarge;
				bt_clen = 2;
				negate = -bdxtail;
				bt_alarge = (double) (negate * ady);
				c = (double) (splitter * negate);
				abig = (double) (c - negate);
				ahi = c - abig;
				alo = negate - ahi;
				c = (double) (splitter * ady);
				abig = (double) (c - ady);
				bhi = c - abig;
				blo = ady - bhi;
				err1 = bt_alarge - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				bt_a[0] = (alo * blo) - err3;
				bt_a[1] = bt_alarge;
				bt_alen = 2;
			} else {
				bdxt_cdy1 = (double) (bdxtail * cdy);
				c = (double) (splitter * bdxtail);
				abig = (double) (c - bdxtail);
				ahi = c - abig;
				alo = bdxtail - ahi;
				c = (double) (splitter * cdy);
				abig = (double) (c - cdy);
				bhi = c - abig;
				blo = cdy - bhi;
				err1 = bdxt_cdy1 - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				bdxt_cdy0 = (alo * blo) - err3;
				bdyt_cdx1 = (double) (bdytail * cdx);
				c = (double) (splitter * bdytail);
				abig = (double) (c - bdytail);
				ahi = c - abig;
				alo = bdytail - ahi;
				c = (double) (splitter * cdx);
				abig = (double) (c - cdx);
				bhi = c - abig;
				blo = cdx - bhi;
				err1 = bdyt_cdx1 - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				bdyt_cdx0 = (alo * blo) - err3;
				_i = (double) (bdxt_cdy0 - bdyt_cdx0);
				bvirt = (double) (bdxt_cdy0 - _i);
				avirt = _i + bvirt;
				bround = bvirt - bdyt_cdx0;
				around = bdxt_cdy0 - avirt;
				bt_c[0] = around + bround;
				_j = (double) (bdxt_cdy1 + _i);
				bvirt = (double) (_j - bdxt_cdy1);
				avirt = _j - bvirt;
				bround = _i - bvirt;
				around = bdxt_cdy1 - avirt;
				_0 = around + bround;
				_i = (double) (_0 - bdyt_cdx1);
				bvirt = (double) (_0 - _i);
				avirt = _i + bvirt;
				bround = bvirt - bdyt_cdx1;
				around = _0 - avirt;
				bt_c[1] = around + bround;
				bt_clarge = (double) (_j + _i);
				bvirt = (double) (bt_clarge - _j);
				avirt = bt_clarge - bvirt;
				bround = _i - bvirt;
				around = _j - avirt;
				bt_c[2] = around + bround;

				bt_c[3] = bt_clarge;
				bt_clen = 4;
				bdyt_adx1 = (double) (bdytail * adx);
				c = (double) (splitter * bdytail);
				abig = (double) (c - bdytail);
				ahi = c - abig;
				alo = bdytail - ahi;
				c = (double) (splitter * adx);
				abig = (double) (c - adx);
				bhi = c - abig;
				blo = adx - bhi;
				err1 = bdyt_adx1 - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				bdyt_adx0 = (alo * blo) - err3;
				bdxt_ady1 = (double) (bdxtail * ady);
				c = (double) (splitter * bdxtail);
				abig = (double) (c - bdxtail);
				ahi = c - abig;
				alo = bdxtail - ahi;
				c = (double) (splitter * ady);
				abig = (double) (c - ady);
				bhi = c - abig;
				blo = ady - bhi;
				err1 = bdxt_ady1 - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				bdxt_ady0 = (alo * blo) - err3;
				_i = (double) (bdyt_adx0 - bdxt_ady0);
				bvirt = (double) (bdyt_adx0 - _i);
				avirt = _i + bvirt;
				bround = bvirt - bdxt_ady0;
				around = bdyt_adx0 - avirt;
				bt_a[0] = around + bround;
				_j = (double) (bdyt_adx1 + _i);
				bvirt = (double) (_j - bdyt_adx1);
				avirt = _j - bvirt;
				bround = _i - bvirt;
				around = bdyt_adx1 - avirt;
				_0 = around + bround;
				_i = (double) (_0 - bdxt_ady1);
				bvirt = (double) (_0 - _i);
				avirt = _i + bvirt;
				bround = bvirt - bdxt_ady1;
				around = _0 - avirt;
				bt_a[1] = around + bround;
				bt_alarge = (double) (_j + _i);
				bvirt = (double) (bt_alarge - _j);
				avirt = bt_alarge - bvirt;
				bround = _i - bvirt;
				around = _j - avirt;
				bt_a[2] = around + bround;

				bt_a[3] = bt_alarge;
				bt_alen = 4;
			}
		}
		if (cdxtail == 0.0) {
			if (cdytail == 0.0) {
				ct_a[0] = 0.0;
				ct_alen = 1;
				ct_b[0] = 0.0;
				ct_blen = 1;
			} else {
				negate = -cdytail;
				ct_alarge = (double) (negate * adx);
				c = (double) (splitter * negate);
				abig = (double) (c - negate);
				ahi = c - abig;
				alo = negate - ahi;
				c = (double) (splitter * adx);
				abig = (double) (c - adx);
				bhi = c - abig;
				blo = adx - bhi;
				err1 = ct_alarge - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				ct_a[0] = (alo * blo) - err3;
				ct_a[1] = ct_alarge;
				ct_alen = 2;
				ct_blarge = (double) (cdytail * bdx);
				c = (double) (splitter * cdytail);
				abig = (double) (c - cdytail);
				ahi = c - abig;
				alo = cdytail - ahi;
				c = (double) (splitter * bdx);
				abig = (double) (c - bdx);
				bhi = c - abig;
				blo = bdx - bhi;
				err1 = ct_blarge - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				ct_b[0] = (alo * blo) - err3;
				ct_b[1] = ct_blarge;
				ct_blen = 2;
			}
		} else {
			if (cdytail == 0.0) {
				ct_alarge = (double) (cdxtail * ady);
				c = (double) (splitter * cdxtail);
				abig = (double) (c - cdxtail);
				ahi = c - abig;
				alo = cdxtail - ahi;
				c = (double) (splitter * ady);
				abig = (double) (c - ady);
				bhi = c - abig;
				blo = ady - bhi;
				err1 = ct_alarge - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				ct_a[0] = (alo * blo) - err3;
				ct_a[1] = ct_alarge;
				ct_alen = 2;
				negate = -cdxtail;
				ct_blarge = (double) (negate * bdy);
				c = (double) (splitter * negate);
				abig = (double) (c - negate);
				ahi = c - abig;
				alo = negate - ahi;
				c = (double) (splitter * bdy);
				abig = (double) (c - bdy);
				bhi = c - abig;
				blo = bdy - bhi;
				err1 = ct_blarge - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				ct_b[0] = (alo * blo) - err3;
				ct_b[1] = ct_blarge;
				ct_blen = 2;
			} else {
				cdxt_ady1 = (double) (cdxtail * ady);
				c = (double) (splitter * cdxtail);
				abig = (double) (c - cdxtail);
				ahi = c - abig;
				alo = cdxtail - ahi;
				c = (double) (splitter * ady);
				abig = (double) (c - ady);
				bhi = c - abig;
				blo = ady - bhi;
				err1 = cdxt_ady1 - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				cdxt_ady0 = (alo * blo) - err3;
				cdyt_adx1 = (double) (cdytail * adx);
				c = (double) (splitter * cdytail);
				abig = (double) (c - cdytail);
				ahi = c - abig;
				alo = cdytail - ahi;
				c = (double) (splitter * adx);
				abig = (double) (c - adx);
				bhi = c - abig;
				blo = adx - bhi;
				err1 = cdyt_adx1 - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				cdyt_adx0 = (alo * blo) - err3;
				_i = (double) (cdxt_ady0 - cdyt_adx0);
				bvirt = (double) (cdxt_ady0 - _i);
				avirt = _i + bvirt;
				bround = bvirt - cdyt_adx0;
				around = cdxt_ady0 - avirt;
				ct_a[0] = around + bround;
				_j = (double) (cdxt_ady1 + _i);
				bvirt = (double) (_j - cdxt_ady1);
				avirt = _j - bvirt;
				bround = _i - bvirt;
				around = cdxt_ady1 - avirt;
				_0 = around + bround;
				_i = (double) (_0 - cdyt_adx1);
				bvirt = (double) (_0 - _i);
				avirt = _i + bvirt;
				bround = bvirt - cdyt_adx1;
				around = _0 - avirt;
				ct_a[1] = around + bround;
				ct_alarge = (double) (_j + _i);
				bvirt = (double) (ct_alarge - _j);
				avirt = ct_alarge - bvirt;
				bround = _i - bvirt;
				around = _j - avirt;
				ct_a[2] = around + bround;

				ct_a[3] = ct_alarge;
				ct_alen = 4;
				cdyt_bdx1 = (double) (cdytail * bdx);
				c = (double) (splitter * cdytail);
				abig = (double) (c - cdytail);
				ahi = c - abig;
				alo = cdytail - ahi;
				c = (double) (splitter * bdx);
				abig = (double) (c - bdx);
				bhi = c - abig;
				blo = bdx - bhi;
				err1 = cdyt_bdx1 - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				cdyt_bdx0 = (alo * blo) - err3;
				cdxt_bdy1 = (double) (cdxtail * bdy);
				c = (double) (splitter * cdxtail);
				abig = (double) (c - cdxtail);
				ahi = c - abig;
				alo = cdxtail - ahi;
				c = (double) (splitter * bdy);
				abig = (double) (c - bdy);
				bhi = c - abig;
				blo = bdy - bhi;
				err1 = cdxt_bdy1 - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				cdxt_bdy0 = (alo * blo) - err3;
				_i = (double) (cdyt_bdx0 - cdxt_bdy0);
				bvirt = (double) (cdyt_bdx0 - _i);
				avirt = _i + bvirt;
				bround = bvirt - cdxt_bdy0;
				around = cdyt_bdx0 - avirt;
				ct_b[0] = around + bround;
				_j = (double) (cdyt_bdx1 + _i);
				bvirt = (double) (_j - cdyt_bdx1);
				avirt = _j - bvirt;
				bround = _i - bvirt;
				around = cdyt_bdx1 - avirt;
				_0 = around + bround;
				_i = (double) (_0 - cdxt_bdy1);
				bvirt = (double) (_0 - _i);
				avirt = _i + bvirt;
				bround = bvirt - cdxt_bdy1;
				around = _0 - avirt;
				ct_b[1] = around + bround;
				ct_blarge = (double) (_j + _i);
				bvirt = (double) (ct_blarge - _j);
				avirt = ct_blarge - bvirt;
				bround = _i - bvirt;
				around = _j - avirt;
				ct_b[2] = around + bround;

				ct_b[3] = ct_blarge;
				ct_blen = 4;
			}
		}

		bctlen = fast_expansion_sum_zeroelim(bt_clen, bt_c, ct_blen, ct_b, bct);
		wlength = scale_expansion_zeroelim(bctlen, bct, adz, w);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother);
		finswap = finnow;
		finnow = finother;
		finother = finswap;

		catlen = fast_expansion_sum_zeroelim(ct_alen, ct_a, at_clen, at_c, cat);
		wlength = scale_expansion_zeroelim(catlen, cat, bdz, w);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother);
		finswap = finnow;
		finnow = finother;
		finother = finswap;

		abtlen = fast_expansion_sum_zeroelim(at_blen, at_b, bt_alen, bt_a, abt);
		wlength = scale_expansion_zeroelim(abtlen, abt, cdz, w);
		finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother);
		finswap = finnow;
		finnow = finother;
		finother = finswap;

		if (adztail != 0.0) {
			vlength = scale_expansion_zeroelim(4, bc, adztail, v);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v, finother);
			finswap = finnow;
			finnow = finother;
			finother = finswap;
		}
		if (bdztail != 0.0) {
			vlength = scale_expansion_zeroelim(4, ca, bdztail, v);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v, finother);
			finswap = finnow;
			finnow = finother;
			finother = finswap;
		}
		if (cdztail != 0.0) {
			vlength = scale_expansion_zeroelim(4, ab, cdztail, v);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v, finother);
			finswap = finnow;
			finnow = finother;
			finother = finswap;
		}

		if (adxtail != 0.0) {
			if (bdytail != 0.0) {
				adxt_bdyt1 = (double) (adxtail * bdytail);
				c = (double) (splitter * adxtail);
				abig = (double) (c - adxtail);
				ahi = c - abig;
				alo = adxtail - ahi;
				c = (double) (splitter * bdytail);
				abig = (double) (c - bdytail);
				bhi = c - abig;
				blo = bdytail - bhi;
				err1 = adxt_bdyt1 - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				adxt_bdyt0 = (alo * blo) - err3;
				c = (double) (splitter * cdz);
				abig = (double) (c - cdz);
				bhi = c - abig;
				blo = cdz - bhi;
				_i = (double) (adxt_bdyt0 * cdz);
				c = (double) (splitter * adxt_bdyt0);
				abig = (double) (c - adxt_bdyt0);
				ahi = c - abig;
				alo = adxt_bdyt0 - ahi;
				err1 = _i - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				u[0] = (alo * blo) - err3;
				_j = (double) (adxt_bdyt1 * cdz);
				c = (double) (splitter * adxt_bdyt1);
				abig = (double) (c - adxt_bdyt1);
				ahi = c - abig;
				alo = adxt_bdyt1 - ahi;
				err1 = _j - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				_0 = (alo * blo) - err3;
				_k = (double) (_i + _0);
				bvirt = (double) (_k - _i);
				avirt = _k - bvirt;
				bround = _0 - bvirt;
				around = _i - avirt;
				u[1] = around + bround;
				u3 = (double) (_j + _k);
				bvirt = u3 - _j;
				u[2] = _k - bvirt;
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
				finswap = finnow;
				finnow = finother;
				finother = finswap;
				if (cdztail != 0.0) {
					c = (double) (splitter * cdztail);
					abig = (double) (c - cdztail);
					bhi = c - abig;
					blo = cdztail - bhi;
					_i = (double) (adxt_bdyt0 * cdztail);
					c = (double) (splitter * adxt_bdyt0);
					abig = (double) (c - adxt_bdyt0);
					ahi = c - abig;
					alo = adxt_bdyt0 - ahi;
					err1 = _i - (ahi * bhi);
					err2 = err1 - (alo * bhi);
					err3 = err2 - (ahi * blo);
					u[0] = (alo * blo) - err3;
					_j = (double) (adxt_bdyt1 * cdztail);
					c = (double) (splitter * adxt_bdyt1);
					abig = (double) (c - adxt_bdyt1);
					ahi = c - abig;
					alo = adxt_bdyt1 - ahi;
					err1 = _j - (ahi * bhi);
					err2 = err1 - (alo * bhi);
					err3 = err2 - (ahi * blo);
					_0 = (alo * blo) - err3;
					_k = (double) (_i + _0);
					bvirt = (double) (_k - _i);
					avirt = _k - bvirt;
					bround = _0 - bvirt;
					around = _i - avirt;
					u[1] = around + bround;
					u3 = (double) (_j + _k);
					bvirt = u3 - _j;
					u[2] = _k - bvirt;
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
					finswap = finnow;
					finnow = finother;
					finother = finswap;
				}
			}
			if (cdytail != 0.0) {
				negate = -adxtail;
				adxt_cdyt1 = (double) (negate * cdytail);
				c = (double) (splitter * negate);
				abig = (double) (c - negate);
				ahi = c - abig;
				alo = negate - ahi;
				c = (double) (splitter * cdytail);
				abig = (double) (c - cdytail);
				bhi = c - abig;
				blo = cdytail - bhi;
				err1 = adxt_cdyt1 - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				adxt_cdyt0 = (alo * blo) - err3;
				c = (double) (splitter * bdz);
				abig = (double) (c - bdz);
				bhi = c - abig;
				blo = bdz - bhi;
				_i = (double) (adxt_cdyt0 * bdz);
				c = (double) (splitter * adxt_cdyt0);
				abig = (double) (c - adxt_cdyt0);
				ahi = c - abig;
				alo = adxt_cdyt0 - ahi;
				err1 = _i - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				u[0] = (alo * blo) - err3;
				_j = (double) (adxt_cdyt1 * bdz);
				c = (double) (splitter * adxt_cdyt1);
				abig = (double) (c - adxt_cdyt1);
				ahi = c - abig;
				alo = adxt_cdyt1 - ahi;
				err1 = _j - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				_0 = (alo * blo) - err3;
				_k = (double) (_i + _0);
				bvirt = (double) (_k - _i);
				avirt = _k - bvirt;
				bround = _0 - bvirt;
				around = _i - avirt;
				u[1] = around + bround;
				u3 = (double) (_j + _k);
				bvirt = u3 - _j;
				u[2] = _k - bvirt;
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
				finswap = finnow;
				finnow = finother;
				finother = finswap;
				if (bdztail != 0.0) {
					c = (double) (splitter * bdztail);
					abig = (double) (c - bdztail);
					bhi = c - abig;
					blo = bdztail - bhi;
					_i = (double) (adxt_cdyt0 * bdztail);
					c = (double) (splitter * adxt_cdyt0);
					abig = (double) (c - adxt_cdyt0);
					ahi = c - abig;
					alo = adxt_cdyt0 - ahi;
					err1 = _i - (ahi * bhi);
					err2 = err1 - (alo * bhi);
					err3 = err2 - (ahi * blo);
					u[0] = (alo * blo) - err3;
					_j = (double) (adxt_cdyt1 * bdztail);
					c = (double) (splitter * adxt_cdyt1);
					abig = (double) (c - adxt_cdyt1);
					ahi = c - abig;
					alo = adxt_cdyt1 - ahi;
					err1 = _j - (ahi * bhi);
					err2 = err1 - (alo * bhi);
					err3 = err2 - (ahi * blo);
					_0 = (alo * blo) - err3;
					_k = (double) (_i + _0);
					bvirt = (double) (_k - _i);
					avirt = _k - bvirt;
					bround = _0 - bvirt;
					around = _i - avirt;
					u[1] = around + bround;
					u3 = (double) (_j + _k);
					bvirt = u3 - _j;
					u[2] = _k - bvirt;
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
					finswap = finnow;
					finnow = finother;
					finother = finswap;
				}
			}
		}
		if (bdxtail != 0.0) {
			if (cdytail != 0.0) {
				bdxt_cdyt1 = (double) (bdxtail * cdytail);
				c = (double) (splitter * bdxtail);
				abig = (double) (c - bdxtail);
				ahi = c - abig;
				alo = bdxtail - ahi;
				c = (double) (splitter * cdytail);
				abig = (double) (c - cdytail);
				bhi = c - abig;
				blo = cdytail - bhi;
				err1 = bdxt_cdyt1 - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				bdxt_cdyt0 = (alo * blo) - err3;
				c = (double) (splitter * adz);
				abig = (double) (c - adz);
				bhi = c - abig;
				blo = adz - bhi;
				_i = (double) (bdxt_cdyt0 * adz);
				c = (double) (splitter * bdxt_cdyt0);
				abig = (double) (c - bdxt_cdyt0);
				ahi = c - abig;
				alo = bdxt_cdyt0 - ahi;
				err1 = _i - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				u[0] = (alo * blo) - err3;
				_j = (double) (bdxt_cdyt1 * adz);
				c = (double) (splitter * bdxt_cdyt1);
				abig = (double) (c - bdxt_cdyt1);
				ahi = c - abig;
				alo = bdxt_cdyt1 - ahi;
				err1 = _j - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				_0 = (alo * blo) - err3;
				_k = (double) (_i + _0);
				bvirt = (double) (_k - _i);
				avirt = _k - bvirt;
				bround = _0 - bvirt;
				around = _i - avirt;
				u[1] = around + bround;
				u3 = (double) (_j + _k);
				bvirt = u3 - _j;
				u[2] = _k - bvirt;
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
				finswap = finnow;
				finnow = finother;
				finother = finswap;
				if (adztail != 0.0) {
					c = (double) (splitter * adztail);
					abig = (double) (c - adztail);
					bhi = c - abig;
					blo = adztail - bhi;
					_i = (double) (bdxt_cdyt0 * adztail);
					c = (double) (splitter * bdxt_cdyt0);
					abig = (double) (c - bdxt_cdyt0);
					ahi = c - abig;
					alo = bdxt_cdyt0 - ahi;
					err1 = _i - (ahi * bhi);
					err2 = err1 - (alo * bhi);
					err3 = err2 - (ahi * blo);
					u[0] = (alo * blo) - err3;
					_j = (double) (bdxt_cdyt1 * adztail);
					c = (double) (splitter * bdxt_cdyt1);
					abig = (double) (c - bdxt_cdyt1);
					ahi = c - abig;
					alo = bdxt_cdyt1 - ahi;
					err1 = _j - (ahi * bhi);
					err2 = err1 - (alo * bhi);
					err3 = err2 - (ahi * blo);
					_0 = (alo * blo) - err3;
					_k = (double) (_i + _0);
					bvirt = (double) (_k - _i);
					avirt = _k - bvirt;
					bround = _0 - bvirt;
					around = _i - avirt;
					u[1] = around + bround;
					u3 = (double) (_j + _k);
					bvirt = u3 - _j;
					u[2] = _k - bvirt;
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
					finswap = finnow;
					finnow = finother;
					finother = finswap;
				}
			}
			if (adytail != 0.0) {
				negate = -bdxtail;
				bdxt_adyt1 = (double) (negate * adytail);
				c = (double) (splitter * negate);
				abig = (double) (c - negate);
				ahi = c - abig;
				alo = negate - ahi;
				c = (double) (splitter * adytail);
				abig = (double) (c - adytail);
				bhi = c - abig;
				blo = adytail - bhi;
				err1 = bdxt_adyt1 - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				bdxt_adyt0 = (alo * blo) - err3;
				c = (double) (splitter * cdz);
				abig = (double) (c - cdz);
				bhi = c - abig;
				blo = cdz - bhi;
				_i = (double) (bdxt_adyt0 * cdz);
				c = (double) (splitter * bdxt_adyt0);
				abig = (double) (c - bdxt_adyt0);
				ahi = c - abig;
				alo = bdxt_adyt0 - ahi;
				err1 = _i - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				u[0] = (alo * blo) - err3;
				_j = (double) (bdxt_adyt1 * cdz);
				c = (double) (splitter * bdxt_adyt1);
				abig = (double) (c - bdxt_adyt1);
				ahi = c - abig;
				alo = bdxt_adyt1 - ahi;
				err1 = _j - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				_0 = (alo * blo) - err3;
				_k = (double) (_i + _0);
				bvirt = (double) (_k - _i);
				avirt = _k - bvirt;
				bround = _0 - bvirt;
				around = _i - avirt;
				u[1] = around + bround;
				u3 = (double) (_j + _k);
				bvirt = u3 - _j;
				u[2] = _k - bvirt;
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
				finswap = finnow;
				finnow = finother;
				finother = finswap;
				if (cdztail != 0.0) {
					c = (double) (splitter * cdztail);
					abig = (double) (c - cdztail);
					bhi = c - abig;
					blo = cdztail - bhi;
					_i = (double) (bdxt_adyt0 * cdztail);
					c = (double) (splitter * bdxt_adyt0);
					abig = (double) (c - bdxt_adyt0);
					ahi = c - abig;
					alo = bdxt_adyt0 - ahi;
					err1 = _i - (ahi * bhi);
					err2 = err1 - (alo * bhi);
					err3 = err2 - (ahi * blo);
					u[0] = (alo * blo) - err3;
					_j = (double) (bdxt_adyt1 * cdztail);
					c = (double) (splitter * bdxt_adyt1);
					abig = (double) (c - bdxt_adyt1);
					ahi = c - abig;
					alo = bdxt_adyt1 - ahi;
					err1 = _j - (ahi * bhi);
					err2 = err1 - (alo * bhi);
					err3 = err2 - (ahi * blo);
					_0 = (alo * blo) - err3;
					_k = (double) (_i + _0);
					bvirt = (double) (_k - _i);
					avirt = _k - bvirt;
					bround = _0 - bvirt;
					around = _i - avirt;
					u[1] = around + bround;
					u3 = (double) (_j + _k);
					bvirt = u3 - _j;
					u[2] = _k - bvirt;
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
					finswap = finnow;
					finnow = finother;
					finother = finswap;
				}
			}
		}
		if (cdxtail != 0.0) {
			if (adytail != 0.0) {
				cdxt_adyt1 = (double) (cdxtail * adytail);
				c = (double) (splitter * cdxtail);
				abig = (double) (c - cdxtail);
				ahi = c - abig;
				alo = cdxtail - ahi;
				c = (double) (splitter * adytail);
				abig = (double) (c - adytail);
				bhi = c - abig;
				blo = adytail - bhi;
				err1 = cdxt_adyt1 - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				cdxt_adyt0 = (alo * blo) - err3;
				c = (double) (splitter * bdz);
				abig = (double) (c - bdz);
				bhi = c - abig;
				blo = bdz - bhi;
				_i = (double) (cdxt_adyt0 * bdz);
				c = (double) (splitter * cdxt_adyt0);
				abig = (double) (c - cdxt_adyt0);
				ahi = c - abig;
				alo = cdxt_adyt0 - ahi;
				err1 = _i - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				u[0] = (alo * blo) - err3;
				_j = (double) (cdxt_adyt1 * bdz);
				c = (double) (splitter * cdxt_adyt1);
				abig = (double) (c - cdxt_adyt1);
				ahi = c - abig;
				alo = cdxt_adyt1 - ahi;
				err1 = _j - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				_0 = (alo * blo) - err3;
				_k = (double) (_i + _0);
				bvirt = (double) (_k - _i);
				avirt = _k - bvirt;
				bround = _0 - bvirt;
				around = _i - avirt;
				u[1] = around + bround;
				u3 = (double) (_j + _k);
				bvirt = u3 - _j;
				u[2] = _k - bvirt;
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
				finswap = finnow;
				finnow = finother;
				finother = finswap;
				if (bdztail != 0.0) {
					c = (double) (splitter * bdztail);
					abig = (double) (c - bdztail);
					bhi = c - abig;
					blo = bdztail - bhi;
					_i = (double) (cdxt_adyt0 * bdztail);
					c = (double) (splitter * cdxt_adyt0);
					abig = (double) (c - cdxt_adyt0);
					ahi = c - abig;
					alo = cdxt_adyt0 - ahi;
					err1 = _i - (ahi * bhi);
					err2 = err1 - (alo * bhi);
					err3 = err2 - (ahi * blo);
					u[0] = (alo * blo) - err3;
					_j = (double) (cdxt_adyt1 * bdztail);
					c = (double) (splitter * cdxt_adyt1);
					abig = (double) (c - cdxt_adyt1);
					ahi = c - abig;
					alo = cdxt_adyt1 - ahi;
					err1 = _j - (ahi * bhi);
					err2 = err1 - (alo * bhi);
					err3 = err2 - (ahi * blo);
					_0 = (alo * blo) - err3;
					_k = (double) (_i + _0);
					bvirt = (double) (_k - _i);
					avirt = _k - bvirt;
					bround = _0 - bvirt;
					around = _i - avirt;
					u[1] = around + bround;
					u3 = (double) (_j + _k);
					bvirt = u3 - _j;
					u[2] = _k - bvirt;
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
					finswap = finnow;
					finnow = finother;
					finother = finswap;
				}
			}
			if (bdytail != 0.0) {
				negate = -cdxtail;
				cdxt_bdyt1 = (double) (negate * bdytail);
				c = (double) (splitter * negate);
				abig = (double) (c - negate);
				ahi = c - abig;
				alo = negate - ahi;
				c = (double) (splitter * bdytail);
				abig = (double) (c - bdytail);
				bhi = c - abig;
				blo = bdytail - bhi;
				err1 = cdxt_bdyt1 - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				cdxt_bdyt0 = (alo * blo) - err3;
				c = (double) (splitter * adz);
				abig = (double) (c - adz);
				bhi = c - abig;
				blo = adz - bhi;
				_i = (double) (cdxt_bdyt0 * adz);
				c = (double) (splitter * cdxt_bdyt0);
				abig = (double) (c - cdxt_bdyt0);
				ahi = c - abig;
				alo = cdxt_bdyt0 - ahi;
				err1 = _i - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				u[0] = (alo * blo) - err3;
				_j = (double) (cdxt_bdyt1 * adz);
				c = (double) (splitter * cdxt_bdyt1);
				abig = (double) (c - cdxt_bdyt1);
				ahi = c - abig;
				alo = cdxt_bdyt1 - ahi;
				err1 = _j - (ahi * bhi);
				err2 = err1 - (alo * bhi);
				err3 = err2 - (ahi * blo);
				_0 = (alo * blo) - err3;
				_k = (double) (_i + _0);
				bvirt = (double) (_k - _i);
				avirt = _k - bvirt;
				bround = _0 - bvirt;
				around = _i - avirt;
				u[1] = around + bround;
				u3 = (double) (_j + _k);
				bvirt = u3 - _j;
				u[2] = _k - bvirt;
				u[3] = u3;
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
				finswap = finnow;
				finnow = finother;
				finother = finswap;
				if (adztail != 0.0) {
					c = (double) (splitter * adztail);
					abig = (double) (c - adztail);
					bhi = c - abig;
					blo = adztail - bhi;
					_i = (double) (cdxt_bdyt0 * adztail);
					c = (double) (splitter * cdxt_bdyt0);
					abig = (double) (c - cdxt_bdyt0);
					ahi = c - abig;
					alo = cdxt_bdyt0 - ahi;
					err1 = _i - (ahi * bhi);
					err2 = err1 - (alo * bhi);
					err3 = err2 - (ahi * blo);
					u[0] = (alo * blo) - err3;
					_j = (double) (cdxt_bdyt1 * adztail);
					c = (double) (splitter * cdxt_bdyt1);
					abig = (double) (c - cdxt_bdyt1);
					ahi = c - abig;
					alo = cdxt_bdyt1 - ahi;
					err1 = _j - (ahi * bhi);
					err2 = err1 - (alo * bhi);
					err3 = err2 - (ahi * blo);
					_0 = (alo * blo) - err3;
					_k = (double) (_i + _0);
					bvirt = (double) (_k - _i);
					avirt = _k - bvirt;
					bround = _0 - bvirt;
					around = _i - avirt;
					u[1] = around + bround;
					u3 = (double) (_j + _k);
					bvirt = u3 - _j;
					u[2] = _k - bvirt;
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u, finother);
					finswap = finnow;
					finnow = finother;
					finother = finswap;
				}
			}
		}

		if (adztail != 0.0) {
			wlength = scale_expansion_zeroelim(bctlen, bct, adztail, w);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother);
			finswap = finnow;
			finnow = finother;
			finother = finswap;
		}
		if (bdztail != 0.0) {
			wlength = scale_expansion_zeroelim(catlen, cat, bdztail, w);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother);
			finswap = finnow;
			finnow = finother;
			finother = finswap;
		}
		if (cdztail != 0.0) {
			wlength = scale_expansion_zeroelim(abtlen, abt, cdztail, w);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w, finother);
			finswap = finnow;
			finnow = finother;
			finother = finswap;
		}

		return finnow[finlength - 1];
	}
}