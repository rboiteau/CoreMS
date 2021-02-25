import numpy as np
import scipy.stats
'''exploratory module based on  Yuanyue Li code at 
TODO add GitHub and Paper here'''

def entropy_distance(v, y):
    merged = v + y
    entropy_increase = 2 * scipy.stats.entropy(merged) - scipy.stats.entropy(v) - scipy.stats.entropy(y)
    return entropy_increase

def _weight_intensity_for_entropy(x):

    if sum(x) > 0:
        WEIGHT_START = 0.25
        WEIGHT_SLOPE = 0.5

        entropy_x = scipy.stats.entropy(x)
        weight = WEIGHT_START + WEIGHT_SLOPE * entropy_x
        x = np.power(x, weight)
        x = x / sum(x)
        return x


def weighted_entropy_distance(v, y):
    v = _weight_intensity_for_entropy(v)
    y = _weight_intensity_for_entropy(y)

    merged = v + y
    entropy_increase = 2 * scipy.stats.entropy(merged) - scipy.stats.entropy(v) - scipy.stats.entropy(y)
    return entropy_increase


def chebyshev_distance(v, y):
    r"""
    Chebyshev distance:

    .. math::

        \underset{i}{\max}{(|v_{i}\ -\ y_{i}|)}
    """
    return np.max(np.abs(v - y))


def squared_euclidean_distance(v, y):
    r"""
    Squared Euclidean distance:

    .. math::

        \sum(v_{i}-y_{i})^2
    """
    return np.sum(np.power(v - y, 2))


def fidelity_similarity(v, y):
    r"""
    Fidelity similarity:

    .. math::

        \sum\sqrt{v_{i}y_{i}}
    """
    return np.sum(np.sqrt(v * y))


def matusita_distance(v, y):
    r"""
    Matusita distance:

    .. math::

        \sqrt{\sum(\sqrt{v_{i}}-\sqrt{y_{i}})^2}
    """
    return np.sqrt(np.sum(np.power(np.sqrt(v) - np.sqrt(y), 2)))


def squared_chord_distance(v, y):
    r"""
    Squared-chord distance:

    .. math::

        \sum(\sqrt{v_{i}}-\sqrt{y_{i}})^2
    """
    return np.sum(np.power(np.sqrt(v) - np.sqrt(y), 2))


def bhattacharya_1_distance(v, y):
    r"""
    Bhattacharya 1 distance:

    .. math::

        (\arccos{(\sum\sqrt{v_{i}y_{i}})})^2
    """
    s = np.sum(np.sqrt(v * y))
    # TODO:Fix this!
    if s > 1:
        if s > 1 + 1e-6:
            print("Error in calculating Bhattacharya 1 distance, got arccos {}".format(s))
        s = 1
    return np.power(np.arccos(s), 2)


def bhattacharya_2_distance(v, y):
    r"""
    Bhattacharya 2 distance:

    .. math::

        -\ln{(\sum\sqrt{v_{i}y_{i}})}
    """
    s = np.sum(np.sqrt(v * y))
    if s == 0:
        return np.inf
    else:
        return -np.log(s)


def harmonic_mean_similarity(v, y):
    r"""
    Harmonic mean similarity:

    .. math::

        #1-2\sum(\frac{v_{i}y_{i}}{v_{i}+y_{i}})
        2\sum(\frac{v_{i}y_{i}}{v_{i}+y_{i}})
    """
    #return 1 - 2 * np.sum(v * y / (v + y))
    return 2 * np.sum(v * y / (v + y))


#def pearson_chi_squared_distance(v, y):
#    r"""
#    Pearson χ2 distance:
#
#    .. math::
#
#        \sum\frac{(v_{i}-y_{i})^2}{y_{i}}
#    """
#    return np.sum(np.power(v - y, 2) / y)


#def neyman_chi_squared_distance(v, y):
#    r"""
#    Neyman χ2 distance:
#
#    .. math::
#
#        \sum\frac{(v_{i}-y_{i})^2}{v_{i}}
#    """
#    return np.sum(np.power(v - y, 2) / v)


#def probabilistic_symmetric_chi_squared_distance(v, y):
#    r"""
#    Probabilistic symmetric χ2 distance:
#
#    .. math::
#
#        \frac{1}{2} \times \sum\frac{(v_{i}-y_{i}\ )^2}{v_{i}+y_{i}\ }
#    """
#    return 1 / 2 * np.sum(np.power(v - y, 2) / (v + y))


#def topsoe_distance(v, y):
#    r"""
#    Topsøe distance:
#
#    .. math::
#
#        \sum{(v_{i}ln\frac{v_{i}}{Z_i}+y_{i}ln\frac{y_{i}}{Z_i}),\ \ \ Z_i=\frac{1}{2}(v_{i}+y_{i})}
#    """
#    z = 1 / 2 * (v + y)
#    z[z == 0] = 1
#    vz = v / z
#    yz = y / z
#    vz[v == 0] = 1
#    yz[y == 0] = 1
#    return np.sum(v * np.log(vz) + y * np.log(yz))


def chernoff_distance(v, y):
    r"""
    Chernoff distance:

    .. math::

        \max{(-ln\sum(v_{i}^ty_{i}^{1-t})^{1-t})},\ t=0.1,\ 0\le\ t<1
    """
    t = 0.1
    return np.max(-np.log(
        np.sum(np.power(np.power(v, t) * np.power(y, 1 - t), 1 - t))))


def ruzicka_distance(v, y):
    r"""
    Ruzicka distance:

    .. math::

        \frac{\sum{|v_{i}-y_{i}|}}{\sum{\max(v_{i},y_{i})}}
    """
    dist = np.sum(np.abs(v - y)) / np.sum(np.maximum(v, y))
    return dist


def roberts_distance(v, y):
    r"""
    Roberts distance:

    .. math::

        1-\sum\frac{(v_{i}+y_{i})\frac{\min{(v_{i},y_{i})}}{\max{(v_{i},y_{i})}}}{\sum(v_{i}+y_{i})}
    """
    return 1 - np.sum((v + y) * np.minimum(v, y) / np.maximum(v, y) / np.sum(v + y))


def intersection_distance(v, y):
    r"""
    Intersection distance:

    .. math::

        1-\frac{\sum\min{(v_{i},y_{i})}}{\min(\sum{v_{i},\sum{y_{i})}}}
    """
    return 1 - np.sum(np.minimum(v, y)) / min(np.sum(v), np.sum(y))


def motyka_distance(v, y):
    r"""
    Motyka distance:

    .. math::

        -\frac{\sum\min{(y_{i},v_{i})}}{\sum(y_{i}+v_{i})}
    """
    dist = np.sum(np.minimum(v, y)) / np.sum(v + y)
    return dist


def canberra_distance(v, y):
    r"""
    Canberra distance:

    .. math::

        #\sum\frac{|v_{i}-y_{i}|}{|v_{i}|+|y_{i}|}
        \sum_{i}\frac{|y_{i} - v_{i}|}{y_{i} + v_{i}}
    """
    #return np.sum(np.abs(v - y) / (np.abs(v) + np.abs(y)))
    return np.sum(np.abs(y - v)/(y + v))

def canberra_metric(v, y):
    r"""
    Canberra Metric

    .. math::

        \frac{1}{\sum_{i}I(v_{i}\neq 0)}\sum_{i}\frac{|y_{i}-v_{i}|}{(y_{i}+v_{i})}
    """

    return (1 / np.sum(v > 0)) * np.sum(np.abs(y - v)/(y + v))


def kulczynski_1_distance(v, y):
    r"""
    Kulczynski 1 distance:

    .. math::

        \frac{\sum{|v_i}-y_i|}{\sum m\ i\ n\ (v_i,y_i)}
    """
    return np.sum(np.abs(y - v)) / np.sum(np.minimum(y, v))


def baroni_urbani_buser_distance(v, y):
    r"""
    Baroni-Urbani-Buser distance:

    .. math::

        1-\frac{\sum\min{(v_i,y_i)}+\sqrt{\sum\min{(v_i,y_i)}\sum(\max{(v)}-\max{(v_i,y_i)})}}{\sum{\max{(v_i,y_i)}+\sqrt{\sum{\min{(v_i,y_i)}\sum(\max{(v)}-\max{(v_i,y_i)})}}}}
    """
    if np.max(v) < np.max(y):
        v, y = y, v
    d1 = np.sqrt(np.sum(np.minimum(v, y) * np.sum(max(v) - np.maximum(v, y))))
    return 1 - (np.sum(np.minimum(v, y)) + d1) / (np.sum(np.maximum(v, y)) + d1)


def penrose_size_distance(v, y):
    r"""
    Penrose size distance:

    .. math::

        \sqrt N\sum{|y_i-v_i|}
    """
    n = np.sum(v > 0)
    return np.sqrt(n) * np.sum(np.abs(y - v))


def mean_character_distance(v, y):
    r"""
    Mean character distance:

    .. math::

        \frac{1}{N}\sum{|y_i-v_i|}
    """
    n = np.sum(v > 0)
    return 1 / n * np.sum(np.abs(y - v))


def lorentzian_distance(v, y):
    r"""
    Lorentzian distance:

    .. math::

        \sum{\ln(1+|v_i-y_i|)}
    """
    return np.sum(np.log(1 + np.abs(y - v)))


def penrose_shape_distance(v, y):
    r"""
    Penrose shape distance:

    .. math::

        \sqrt{\sum((v_i-\bar{v})-(y_i-\bar{y}))^2}
    """
    v_avg = np.mean(v)
    y_avg = np.mean(y)
    return np.sqrt(np.sum(np.power((y - y_avg) - (v - v_avg), 2)))


def clark_distance(v, y):
    r"""
    Clark distance:

    .. math::

        #(\frac{1}{N}\sum(\frac{v_i-y_i}{|v_i|+|y_i|})^2)^\frac{1}{2}
        \sqrt{\sum(\frac{|v_i-y_i|}{v_i+y_i})^2}
    """
    #n = np.sum(v > 0)
    #return np.sqrt(1 / n * np.sum(np.power((v - y) / (np.abs(v) + np.abs(y)), 2)))
    return np.sqrt(np.sum(np.power(np.abs(y - v) / (y + v), 2)))


def hellinger_distance(v, y):
    r"""
    Hellinger distance:

    .. math::

        #\sqrt{2\sum(\sqrt{\frac{v_i}{\bar{v}}}-\sqrt{\frac{y_i}{\bar{y}}})^2}
        \sqrt{2\sum(\sqrt{v_i}-\sqrt{y_i})^2}
    """
    #v_avg = np.mean(v)
    #y_avg = np.mean(y)
    #return np.sqrt(2 * np.sum(np.power(np.sqrt(v / v_avg) - np.sqrt(y / y_avg), 2)))
    return np.sqrt(2 * np.sum(np.power(np.sqrt(y) - np.sqrt(v), 2)))


def whittaker_index_of_association_distance(v, y):
    r"""
    Whittaker index of association distance:

    .. math::

        \frac{1}{2}\sum|\frac{v_i}{\bar{v}}-\frac{y_i}{\bar{y}}|
    """
    v_avg = np.mean(v)
    y_avg = np.mean(y)
    return 1 / 2 * np.sum(np.abs(v / v_avg - y / y_avg))


#def symmetric_chi_squared_distance(v, y):
#    r"""
#    Symmetric χ2 distance:
#
#    .. math::
#
#        \sqrt{\sum{\frac{\bar{v}+\bar{y}}{N(\bar{v}+\bar{y})^2}\frac{(v_i\bar{y}-y_i\bar{v})^2}{v_i+y_i}\ }}
#    """
#    v_avg = np.mean(v)
#    y_avg = np.mean(y)
#    n = np.sum(v > 0)
#
#    d1 = (v_avg + y_avg) / (n * np.power(v_avg + y_avg, 2))
#    return np.sqrt(d1 * np.sum(np.power(v * y_avg - y * v_avg, 2) / (v + y)))

def similarity_index_distance(v, y):
    r"""
    Similarity Index Distance:

    .. math::

        \sqrt{\frac{\sum\{\frac{v_i-y_i}{y_i}\}^2}{N}}
    """
    n = np.sum(v > 0)
    return np.sqrt(1 / n * np.sum(np.power((v - y) / y, 2)))


def improved_similarity_distance(v, y):
    r"""
    Improved Similarity Index:

    .. math::

        \sqrt{\frac{1}{N}\sum\{\frac{y_i-v_i}{y_i+v_i}\}^2}
    """
    n = np.sum(v > 0)
    return np.sqrt(1 / n * np.sum(np.power((y - v) / (y + v), 2)))


def absolute_value_distance(v, y):
    r"""
    Absolute Value Distance:

    .. math::

        \frac { \sum(|y_i-v_i|)}{\sum v_i}

    """
    dist = np.sum(np.abs(y - v)) / np.sum(v)
    return dist

def spectral_contrast_angle_distance(v, y):
    r"""
    Spectral Contrast Angle:

    .. math::

        1 - \frac{\sum{y_iv_i}}{\sqrt{\sum y_i^2\sum v_i^2}}
        \arccos(\frac{\sum_{P}y_{p}^* v_{p}^*}{\sqrt{\sum_{P}y_{p}^{*2} \sum_{P}v_{p}^{*2}}})
    """
    #return 1 - np.sum(y * v) / \
    #       np.sqrt(np.sum(np.power(y, 2)) * np.sum(np.power(v, 2)))

    return np.arccos(np.sum(y * v) / (np.sqrt(np.sum(np.power(y, 2)) * np.sum(np.power(v, 2)))))


def wave_hedges_distance(v, y):
    r"""
    Wave Hedges distance:

    .. math::

        \sum\frac{|v_i-y_i|}{\max{(v_i,y_i)}}
    """
    return np.sum(np.abs(v - y) / np.maximum(v, y))

def dice_similarity(v, y):
    r"""
    Dice similarity:

    .. math::

        \frac{\sum(v_i-y_i)^2}{\sum v_i^2+\sum y_i^2}
        \frac{2 * \sum_{i}v_{i}y_{i}}{\sum_{i}y_{i}^2 + \sum_{i}v_{i}^2}
    """
    return 2 * np.sum(v * y) / (np.sum(np.power(v, 2)) + np.sum(np.power(y, 2)))


def inner_product_distance(v, y):
    r"""
    Inner Product distance:

    .. math::

        1-\sum{v_iy_i}
    """
    return 1 - np.sum(v * y)


def divergence_distance(v, y):
    r"""
    Divergence distance:

    .. math::

        2\sum\frac{(v_i-y_i)^2}{(v_i+y_i)^2}
    """
    return 2 * np.sum((np.power(v - y, 2)) / np.power(v + y, 2))


def _chi_squared_distance(v, y):
    r"""
    Additive symmetric χ2 distance:

    .. math::

        \sum\frac{(v_i-y_i)^2(v_i+y_i)}{v_iy_i}
    """
    dist = np.sum(np.power(v - y, 2) * (v + y) / (v * y))
    return dist


def jensen_difference_distance(v, y):
    r"""
    Jensen difference:

    .. math::

        \sum[\frac{1}{2}(v_i\ln{v_i}+y_i\ln{y_i})-(\frac{v_i+y_i}{2})\ln{(\frac{v_i+y_i}{2})}]
    """
    y_v_avg = (y + v) / 2
    return np.sum(
        1 / 2 * (y * np.log(y) + v * np.log(v)) -
        y_v_avg * np.log(y_v_avg)
    )


def kumar_johnson_distance(v, y):
    r"""
    Kumar-Johnson distance:

    .. math::

        \sum\frac{(v_i^2-y_i^2)^2}{2(v_iy_i)^\frac{3}{2}}
    """
    return np.sum(
        np.power(np.power(v, 2) - np.power(y, 2), 2) / \
        (2 * np.power(v * y, 3 / 2))
    )


def avg_l_distance(v, y):
    r"""
    Avg (L1, L∞) distance:

    .. math::

        \frac{1}{2}(\sum|v_i-y_i|+\underset{i}{\max}{|v_i-y_i|})
    """
    return 1 / 2 * (np.sum(np.abs(v - y)) + max(np.abs(v - y)))


def vicis_wave_hadges_distance(v, y):
    r"""
    Vicis-Wave Hadges distance:

    .. math::

        \sum\frac{|v_i-y_i|}{\min{(v_i,\ y_i)}}
    """
    return np.sum(np.abs(v - y) / np.minimum(v, y))


def vicis_symmetric_chi_squared_1_distance(v, y):
    r"""
    Vicis-Symmetric χ2 1 distance:

    .. math::

        \sum\frac{(v_i-y_i)^2}{\min{(v_i,y_i)^2}}
    """
    return np.sum(np.power(v - y, 2) / np.power(np.minimum(v, y), 2))


def vicis_symmetric_chi_squared_2_distance(v, y):
    r"""
    Vicis-Symmetric χ2 2 distance:

    .. math::

        \sum\frac{(v_i-y_i)^2}{\min{(v_i,y_i)}}
    """
    return np.sum(np.power(v - y, 2) / np.minimum(v, y))


def vicis_symmetric_chi_squared_3_distance(v, y):
    r"""
    Vicis-Symmetric χ2 3 distance:

    .. math::

        \sum\frac{(v_i-y_i)^2}{\max{(v_i,y_i)}}
    """
    return np.sum(np.power(v - y, 2) / np.maximum(v, y))


def max_symmetric_chi_squared_distance(v, y):
    r"""
    Max-Symmetric χ2 distance:

    .. math::

        \max{(\sum\frac{(v_i-y_i)^2}{v_i},\sum\frac{(v_i-y_i)^2}{y_i})}
    """
    return max(np.sum(np.power(v - y, 2) / v), np.sum(np.power(v - y, 2) / y))


def min_symmetric_chi_squared_distance(v, y):
    r"""
    Min-Symmetric χ2 distance:

    .. math::

        \min{(\sum\frac{(v_i-y_i)^2}{v_i},\sum\frac{(v_i-y_i)^2}{y_i})}
    """
    return min(np.sum(np.power(v - y, 2) / v), np.sum(np.power(v - y, 2) / y))

"""added by Allison"""
def additive_sym_chi_sq(v, y):
    r"""
    Additive Symmetric χ2 distance:

    .. math::

        \sum_{i}\frac{(y_{i} - v_{i})^2(y_{i}+v_{i})}{y_{i}v_{i}}
    """
    return np.sum((np.power(y - v, 2) * (y + v))/(y * v))

def bhattacharya_distance(v, y):
    r"""
    Bhattacharya Distance:

    .. math::

        -ln(\sum_{i}\sqrt{y_{i}v_{i}})
    """
    return -1 * np.log(np.sum(np.sqrt(y * v)))

def generalized_ochai_index(v, y):
    r"""
    Generalized Ochai Index

    .. math::

        1 - \frac{\sum_{i}min(y_{i}, v_{i})}{\sqrt{\sum_{i}y_{i} \sum_{i}v_{i}}}
    """

    ind = np.sum(np.minimum(y, v)) / np.sqrt(np.sum(y) * np.sum(v))
    return 1 - ind 

def gower_distance(v, y):
    r"""
    Gower Distance

    .. math::

        \frac{1}{N}\sum_{i}|y_{i} - v_{i}|
    """

    n = np.sum(y > 0)
    return (1 / n) * np.sum(np.abs(y - v))

def impr_sqrt_cosine_sim(v, y):
    r"""
    Improved Square Root Cosine Similarity

    .. math::

        \frac{\sum_{i}\sqrt{y_{i}v_{i}}}{\sum_{i}\sqrt{y_{i}}\sum_{i}\sqrt{v_{i}}}
    """

    return np.sum(np.sqrt(y * v)) / (np.sum(np.sqrt(y)) * np.sum(np.sqrt(v)))

def intersection_sim(v, y):
    r"""
    Intersection Similarity

    .. math::

        \sum_{i}min(y_{i}, v_{i})
    """

    return np.sum(np.minimum(y, v))

def j_divergence(v, y):
    r"""
    J Divergence

    .. math::
        
        \sum_{i}(y_{i} - v_{i}) ln(\frac{y_{i}}{v_{i}})
    """

    return np.sum((v - y) * np.log(v / y))

def jensen_shannon_index(v, y):
    r"""
    Jensen-Shannon Index

    .. math::

        \frac{1}{2}[\sum_{i}y_{i}ln(\frac{2y_{i}}{y_{i} + v_{i}}) + \sum_{i}v_{i}ln(\frac{2v_{i}}{y_{i}+v_{i}})]
    """

    return (1 / 2) * (np.sum(y * np.log(2 * y / (y + v))) + np.sum(v * np.log(2 * v / (y + v))))

def k_divergence(v, y):
    r"""
    K-Divergence

    .. math::

        \sum_{i}y_{i}ln(\frac{2y_{i}}{y_{i} + v_{i}})
    """

    return np.sum(v * np.log((2 * v) / (y + v)))

"""added by Chae"""
def topsoe_distance(v, y):
    r""" Fixed
    "I commented out the previous one; please review"
    """
    #v[v==0] = 1 #added by amt
    #y[y==0] = 1 #added by amt
    return np.sum((y * np.log((2 * y)/(y + v))) + (v * np.log((2 * v)/(y + v))))

def probabilistic_symmetric_chi_squared_distance(v, y):
    r""" Fixed
    "I commented out the previous one; please review"
    """
    return 2 * np.sum(np.sum(np.power(y - v, 2) / (y + v)))

def VW6(v, y):
    r"""
    "appears to be the same as max_symmetric_chi_squared_distance"
    """
    return min(np.sum(np.power(y - v, 2) / y), np.sum(np.power(y - v, 2) / v))

def VW5(v, y):
    r"""
    "appears to be the same as max_symmetric_chi_squared_distance"
    """
    return max(np.sum(np.power(y - v, 2) / y), np.sum(np.power(y - v, 2) / v))

def VW4(v, y):
    r"""
    "Tecnically the Symmetric chi2 eq63"
    """
    return np.sum(np.power(y - v, 2) / np.maximum(y, v))

def VW3(v, y):
    r"""
    "New"
    """
    return np.sum(np.power(y - v, 2) / np.minimum(y, v))

def VW2(v, y):
    r"""
    "New"
    """
    return np.sum(np.power(y - v, 2) / np.power(np.minimum(y, v), 2))

def VW1(v, y):
    r"""
    "New"
    """
    return np.sum(np.abs(y - v) / np.minimum(y, v))

def taneja_divergence(v, y):
    r"""
    "New"
    """
    return np.sum(((y + v) / 2) * np.log((y + v)/(2 * np.sqrt(y * v))))

def symmetric_chi_squared_distance (v, y):
    r"""
    "New"
    """
    return np.sum(np.power(y - v, 2) / (y * v))

def squared_chi_squared_distance(v, y):
    r"""
    "New"
    """
    return np.sum(np.power(y - v, 2) / (y + v))

def square_root_cosine_correlation(v, y):
    r"""
    "New"
    """
    return np.sum(np.sqrt(y * v)) / (np.sum(y) * np.sum(v))

def sorensen_distance(v, y):
    r"""
    "New"
    """
    return np.sum(np.abs(y - v)) / (np.sum(y + v))

def Pearson_chi_squared_distance(v, y):
    r"""
    "New"
    """
    return np.sum(np.power(y - v, 2) / v)

def Neyman_chi_squared_distance(v, y):
    r"""
    "New"
    """
    return np.sum(np.power(y - v, 2) / y)

def Minokowski_3(v, y):
    r"""
    "New"
    """
    return np.power(np.sum(np.power(np.abs(y - v), 3)), 1/3)

def Minokowski_4(v, y):
    r"""
    "New"
    """
    return np.power(np.sum(np.power(np.abs(y - v), 4)), 1/4)

def kumarjohnson_divergence(v, y):
    r"""
    "New"
    """
    return np.sum(np.power(np.power(y, 2) + np.power(v, 2), 2) / (2* np.power(y * v, 3/2)))

def kumarhassebrook_similarity(v, y):
    r"""
    "New"
    """
    return np.sum(y * v) / (np.sum(np.power(y, 2)) + np.sum(np.power(v, 2)) - np.sum(y * v))

def kullbackleibler_divergence (v, y):
    r"""
    "New"
    """
    return np.sum(v * np.log(v / y))

def soergel_distance(v, y):
    r"""
    "New"
    """
    return np.sum(np.abs(y - v))/np.sum(np.maximum(y, v))
