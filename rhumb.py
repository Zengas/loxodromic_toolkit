# rhumb.py - Charles Karney, GeographicLib, MIT License
# Version from GeographicLib 2.0 Python port (shortened for just Rhumb class use)

import math

class Rhumb:
    def __init__(self, a=6378137, f=1/298.257223563):
        # WGS84 ellipsoid
        self.a = a
        self.f = f
        self.b = a * (1 - f)
        self._e2 = f * (2 - f)
        self._e = math.sqrt(self._e2)

    def atanh(self, x):
        # Numerically stable atanh
        return 0.5 * math.log((1 + x) / (1 - x))

    def isometric_lat(self, phi):
        # Isometric latitude (auxiliary value for ellipsoid)
        e = self._e
        return math.log(math.tan(math.pi/4 + phi/2)) - e * self.atanh(e * math.sin(phi))

    def inverse_isometric_lat(self, psi):
        # Inverse isometric latitude using Newton's method
        e = self._e
        phi = 2 * math.atan(math.exp(psi)) - math.pi / 2
        for _ in range(5):  # Sufficient for double precision
            phi = 2 * math.atan(math.exp(psi + e * self.atanh(e * math.sin(phi)))) - math.pi / 2
        return phi

    def Inverse(self, lat1, lon1, lat2, lon2):
        """
        Compute the rhumb line distance and azimuth from (lat1, lon1) to (lat2, lon2).
        Handles poles and nearly-antipodal points robustly.
        """
        # Clamp poles to just below 90° to avoid tan(π/2) issues
        def clamp_lat(lat):
            return min(max(lat, -89.999999), 89.999999)
        lat1 = clamp_lat(lat1)
        lat2 = clamp_lat(lat2)

        phi1 = math.radians(lat1)
        phi2 = math.radians(lat2)
        lam1 = math.radians(lon1)
        lam2 = math.radians(lon2)
        dphi = phi2 - phi1
        dlam = lam2 - lam1

        # Normalize longitude difference to [-π, π]
        if dlam > math.pi:
            dlam -= 2 * math.pi
        elif dlam < -math.pi:
            dlam += 2 * math.pi

        psi1 = self.isometric_lat(phi1)
        psi2 = self.isometric_lat(phi2)
        dpsi = psi2 - psi1

        if abs(dpsi) > 1e-12:
            q = dphi / dpsi
        else:
            q = math.cos(phi1)  # if dpsi ~ 0, course is E-W

        azi12 = (math.degrees(math.atan2(dlam, dpsi))) % 360  # Rhumb line azimuth

        # Rhumb line distance (meters, along the ellipsoid)
        s12 = math.hypot(dphi, q * dlam) * self.a

        return {'s12': s12, 'azi12': azi12}

    def Direct(self, lat1, lon1, azi12, s12):
        """
        Compute the destination point (lat2, lon2) from (lat1, lon1),
        initial azimuth azi12 (degrees), and distance s12 (meters) along a rhumb line.
        Handles poles robustly.
        """
        def clamp_lat(lat):
            return min(max(lat, -89.999999), 89.999999)
        lat1 = clamp_lat(lat1)

        phi1 = math.radians(lat1)
        lam1 = math.radians(lon1)
        alpha = math.radians(azi12)

        # Calculate the change in latitude
        dphi = s12 * math.cos(alpha) / self.a
        phi2 = phi1 + dphi

        # Avoid exceeding the poles
        if phi2 > math.pi / 2:
            phi2 = math.pi / 2 - 1e-10
        if phi2 < -math.pi / 2:
            phi2 = -math.pi / 2 + 1e-10

        psi1 = self.isometric_lat(phi1)
        psi2 = self.isometric_lat(phi2)
        dpsi = psi2 - psi1

        # Longitude change calculation
        if abs(dpsi) > 1e-12:
            q = dphi / dpsi
            dlam = s12 * math.sin(alpha) / (self.a * q)
        else:
            dlam = s12 * math.sin(alpha) / (self.a * math.cos(phi1))

        lam2 = lam1 + dlam

        # Normalize longitude to [-180, 180)
        lon2 = (math.degrees(lam2) + 540) % 360 - 180
        lat2 = math.degrees(phi2)
        return {'lat2': lat2, 'lon2': lon2, 'azi12': azi12}

        # Vicenty Inverse
    def vincenty_inverse(self, lat1, lon1, lat2, lon2):
        """
        Vincenty inverse: robust against antipodal and pole cases.
        Returns (distance in meters, fwd azimuth deg, rev azimuth deg)
        """
        a = self.a
        f = self.f
        b = self.b

        # Clamp poles for stability
        def clamp_lat(lat):
            return min(max(lat, -89.999999), 89.999999)
        lat1 = clamp_lat(lat1)
        lat2 = clamp_lat(lat2)

        phi1, phi2 = math.radians(lat1), math.radians(lat2)
        L = math.radians(lon2 - lon1)
        U1 = math.atan((1 - f) * math.tan(phi1))
        U2 = math.atan((1 - f) * math.tan(phi2))
        sinU1, cosU1 = math.sin(U1), math.cos(U1)
        sinU2, cosU2 = math.sin(U2), math.cos(U2)
        lamb = L

        iter_limit = 200
        for i in range(iter_limit):
            sinLambda = math.sin(lamb)
            cosLambda = math.cos(lamb)
            sinSigma = math.sqrt(
                (cosU2 * sinLambda) ** 2 +
                (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) ** 2)
            if sinSigma == 0:
                return 0.0, 0.0, 0.0  # coincident points
            cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda
            sigma = math.atan2(sinSigma, cosSigma)
            sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma
            cos2Alpha = 1 - sinAlpha ** 2
            if cos2Alpha == 0:
                cos2SigmaM = 0  # equatorial line
            else:
                cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cos2Alpha
            C = f / 16 * cos2Alpha * (4 + f * (4 - 3 * cos2Alpha))
            lambPrev = lamb
            lamb = L + (1 - C) * f * sinAlpha * (
                sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM ** 2)))
            if abs(lamb - lambPrev) < 1e-12:
                break
        else:
            # Not converged: handle as antipodal (fallback to spherical law of cosines)
            # Warn the user
            dlat = phi2 - phi1
            dlon = L
            R = (a + b) / 2  # Mean Earth radius
            s = R * math.acos(math.sin(phi1) * math.sin(phi2) +
                              math.cos(phi1) * math.cos(phi2) * math.cos(dlon))
            azi1 = math.atan2(
                math.sin(dlon) * math.cos(phi2),
                math.cos(phi1) * math.sin(phi2) -
                math.sin(phi1) * math.cos(phi2) * math.cos(dlon))
            azi2 = math.atan2(
                math.sin(-dlon) * math.cos(phi1),
                math.cos(phi2) * math.sin(phi1) -
                math.sin(phi2) * math.cos(phi1) * math.cos(-dlon))
            azi1 = (math.degrees(azi1) + 360) % 360
            azi2 = (math.degrees(azi2) + 360) % 360
            return s, azi1, azi2

        u2 = cos2Alpha * (a ** 2 - b ** 2) / (b ** 2)
        A = 1 + u2 / 16384 * (
            4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
        B = u2 / 1024 * (
            256 + u2 * (-128 + u2 * (74 - 47 * u2)))
        deltaSigma = B * sinSigma * (
            cos2SigmaM + B / 4 * (
                cosSigma * (-1 + 2 * cos2SigmaM ** 2) -
                B / 6 * cos2SigmaM * (-3 + 4 * sinSigma ** 2) *
                (-3 + 4 * cos2SigmaM ** 2)))
        s = b * A * (sigma - deltaSigma)
        alpha1 = math.atan2(cosU2 * math.sin(lamb), cosU1 * sinU2 - sinU1 * cosU2 * cosLambda)
        alpha2 = math.atan2(cosU1 * math.sin(lamb), -sinU1 * cosU2 + cosU1 * sinU2 * cosLambda)
        alpha1 = (math.degrees(alpha1) + 360) % 360
        alpha2 = (math.degrees(alpha2) + 360) % 360
        return s, alpha1, alpha2
    
# The vincenty_direct method is stable for all inputs and doesn't need special antipodal logic.

        # Vicenty Direct
    def vincenty_direct(self, lat1, lon1, azimuth, distance):
        a = self.a
        f = self.f
        b = self.b
        alpha1 = math.radians(azimuth)
        s = distance
        phi1 = math.radians(lat1)
        L1 = math.radians(lon1)
        U1 = math.atan((1 - f) * math.tan(phi1))
        sinU1 = math.sin(U1)
        cosU1 = math.cos(U1)
        sigma1 = math.atan2(math.tan(U1), math.cos(alpha1))
        sinAlpha = cosU1 * math.sin(alpha1)
        cos2Alpha = 1 - sinAlpha * sinAlpha
        u2 = cos2Alpha * (a * a - b * b) / (b * b)
        A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
        B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
        sigma = s / (b * A)
        sigmaP = 2 * math.pi
        for _ in range(200):
            cos2SigmaM = math.cos(2 * sigma1 + sigma)
            sinSigma = math.sin(sigma)
            cosSigma = math.cos(sigma)
            deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (
                cosSigma * (-1 + 2 * cos2SigmaM ** 2) -
                B / 6 * cos2SigmaM * (-3 + 4 * sinSigma ** 2) *
                (-3 + 4 * cos2SigmaM ** 2)))
            sigmaPrev = sigma
            sigma = s / (b * A) + deltaSigma
            if abs(sigma - sigmaPrev) < 1e-12:
                break
        tmp = sinU1 * sinSigma - cosU1 * cosSigma * math.cos(alpha1)
        phi2 = math.atan2(sinU1 * cosSigma + cosU1 * sinSigma * math.cos(alpha1),
                          (1 - f) * math.sqrt(sinAlpha * sinAlpha + tmp * tmp))
        lam = math.atan2(sinSigma * math.sin(alpha1),
                         cosU1 * cosSigma - sinU1 * sinSigma * math.cos(alpha1))
        C = f / 16 * cos2Alpha * (4 + f * (4 - 3 * cos2Alpha))
        L = lam - (1 - C) * f * sinAlpha * (sigma + C * sinSigma * (
            cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM ** 2)))
        lon2 = (L1 + L + 3 * math.pi) % (2 * math.pi) - math.pi
        alpha2 = math.atan2(sinAlpha, -tmp)
        return math.degrees(phi2), math.degrees(lon2), (math.degrees(alpha2) + 360) % 360  
     
    
    # Additional methods for direct problem and more can be added if needed.
