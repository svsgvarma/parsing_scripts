import os.path
import sys
from optparse import OptionParser
from vcfReader import *
#import pylab
import operator
from itertools import *
import math
import random
    
DEBUG = False

class Range:
    ANY = '*'
    def __init__(self, left = ANY, right = ANY, leftOpen = False, rightOpen = False):
        self.left = left
        self.right = right
        self.leftOpen = leftOpen
        self.rightOpen = rightOpen
        
    def __str__(self):
        leftB, rightB = '[', ']'
        if self.leftOpen: leftB = '('
        if self.rightOpen: rightB = ')'
        return '%s%s, %s%s' % (leftB, self.left, self.right, rightB)
        
    __repr__ = __str__
        
    def dashedString(self):
        return str(self).replace(", ", "-")
        
    def contains(self, v):
        def test(r, op, open):
            return r == Range.ANY or op(v, r) or (not open and v == r)  
        return test(self.left, operator.__gt__, self.leftOpen) and test(self.right, operator.__lt__, self.rightOpen)

class CallCovariate:
    def __init__(self, feature, featureRange, qualRange, FPRate = None):
        self.feature = feature
        
        self.featureRange = featureRange

        self.qualRange = qualRange
        self.FPRate = FPRate
        
    def containsVariant(self, call):
        inFeature = self.featureRange.contains(call.getField(self.feature))
        inQual = self.qualRange.contains(call.getQual())
        #print 'inFeature, inQual', inFeature, inQual
        return inFeature and inQual

    def getFPRate(self): return self.FPRate
    def getFeature(self): return self.feature
    
    def getCovariateField(self): return self.getFeature() + '_RQ'
    
    def __str__(self): return "[CC feature=%s range=%s qualRange=%s]" % (self.feature, self.featureRange, self.qualRange)

class RecalibratedCall:
    def __init__(self, call, features):
        self.call = call
        self.features = dict([[feature, None] for feature in features])
        
    def recalFeature( self, feature, FPRate ):
        assert self.features[feature] == None, "Feature " + feature + ' has value ' + str(self.features[feature]) + ' for call ' + str(self.call) # not reassigning values
        assert FPRate <= 1 and FPRate >= 0
        self.features[feature] = FPRate
        
    def getFeature( self, feature, missingValue = None, phredScaleValue = False ):
        v = self.features[feature]
        if v == None:
            return missingValue
        elif phredScaleValue:
            return phredScale(v)
        else:
            return v
        
    def jointFPErrorRate(self):
        #print self.features
        logTPRates = [math.log10(1-r) for r in self.features.itervalues() if r <> None]
        logJointTPRate = reduce(lambda x, y: x + y, logTPRates, 0)
        logJointTPRate = min(logJointTPRate, 1e-3 / 3) # approximation from het of 0.001
        jointTPRate = math.pow(10, logJointTPRate)
        #print logTPRates
        #print logJointTPRate, jointTPRate
        return 1 - jointTPRate
        
    def featureStringList(self):
        return ','.join(map(lambda feature: '%s=Q%d' % (feature, self.getFeature(feature, '*', True)), self.features.iterkeys()))        
        
    def __str__(self):
        return '[%s: %s => Q%d]' % (str(self.call), self.featureStringList(), phredScale(self.jointFPErrorRate()))

def readVariants( file, maxRecords = None, decodeAll = True, downsampleFraction = 1, filter = None, minQScore = -1, maxQScore = 10000000, mustBeVariant = False, skip = None ):
    if filter == None:
        filter = not OPTIONS.unfiltered
        
    f = open(file)
    header, columnNames, lines = readVCFHeader(f)

    nLowQual = 0
    counter = 0
    def parseVariant(args):
        global nLowQual, counter
        header1, VCF, counter = args
        counter += 1
        #print counter, skip, counter % skip
        if filter and not VCF.passesFilters() or ( False and mustBeVariant == True and not VCF.isVariant() ): # currently ignore mustBeVariant
            #print 'filtering', VCF
            return None
        elif VCF.getQual() <= minQScore or VCF.getQual() > maxQScore:
            #print 'filtering', VCF.getQual()
            #nLowQual += 1
            return None
        elif skip <> None and counter % skip <> 0:
            #print 'skipping'
            return None
        elif random.random() <= downsampleFraction:
            #print 'keeping', VCF
            return VCF
        else:
            return None

    variants = ifilter(None, imap(parseVariant, islice(lines2VCF(lines, header=header, columnNames = columnNames, extendedOutput = True, decodeAll = decodeAll), maxRecords)))
    if nLowQual > 0:
        print '%d snps filtered due to QUAL < %d' % (nLowQual, minQScore)
    return header, variants

def selectVariants( variants, selector = None ):
    if selector <> None:
        return filter(selector, variants)
    else:
        return variants

def titv(variants):
    ti = len(filter(VCFRecord.isTransition, variants))
    tv = len(variants) - ti
    titv = ti / (1.0*max(tv,1))

    return titv

def dbSNPRate(variants):
    inDBSNP = len(filter(VCFRecord.isKnown, variants))
    return float(inDBSNP) / max(len(variants),1)

def gaussian(x, mu, sigma):    
    constant = 1 / math.sqrt(2 * math.pi * sigma**2)
    exponent = -1 * ( x - mu )**2 / (2 * sigma**2)
    return constant * math.exp(exponent)

# if target = T, and FP calls have ti/tv = 0.5, we want to know how many FP calls
# there are in N calls with ti/tv of X.  
# 
def titvFPRateEstimate(variants, target):
    titvRatio = titv(variants)
    
    # f <- function(To,T) { (To - T) / (1/2 - T) + 0.001 }
    def theoreticalCalc():
        if titvRatio >= target:
            FPRate = 0
        else:
            FPRate = (titvRatio - target) / (0.5 - target)
        FPRate = min(max(FPRate, 0), 1)
        TPRate = max(min(1 - FPRate, 1 - dephredScale(OPTIONS.maxQScoreForCovariate)), dephredScale(OPTIONS.maxQScoreForCovariate))
        if DEBUG: print 'FPRate', FPRate, titvRatio, target
        assert FPRate >= 0 and FPRate <= 1
        return TPRate
    
    # gaussian model
    def gaussianModel():
        LEFT_HANDED = True
        sigma = 1 # old value is 5
        constant = 1 / math.sqrt(2 * math.pi * sigma**2)
        exponent = -1 * ( titvRatio - target )**2 / (2 * sigma**2)
        TPRate = gaussian(titvRatio, target, sigma) / gaussian(target, target, sigma)
        if LEFT_HANDED and titvRatio >= target:
            TPRate = 1
        TPRate -= dephredScale(OPTIONS.maxQScoreForCovariate)
        if DEBUG: print 'TPRate', TPRate, constant, exponent, dephredScale(OPTIONS.maxQScoreForCovariate)
        return TPRate
    
    FPRate = 1 - theoreticalCalc()
    #FPRate = 1 - gaussianModel()
    nVariants = len(variants)
    
    if DEBUG: print ':::', nVariants, titvRatio, target, FPRate
    
    return titvRatio, FPRate
    
def phredScale(errorRate):
    return -10 * math.log10(max(errorRate, 1e-10))

def dephredScale(qscore):
    return math.pow(10, float(qscore) / -10)

def frange6(*args):
    """A float range generator."""
    start = 0.0
    step = 1.0

    l = len(args)
    if l == 1:
        end = args[0]
    elif l == 2:
        start, end = args
    elif l == 3:
        start, end, step = args
        if step == 0.0:
            raise ValueError, "step must not be zero"
    else:
        raise TypeError, "frange expects 1-3 arguments, got %d" % l

    v = start
    while True:
        if (step > 0 and v >= end) or (step < 0 and v <= end):
            raise StopIteration
        yield v
        v += step

def compareFieldValues( v1, v2 ):
    if type(v1) <> type(v2):
        #print 'Different types', type(v1), type(v2)
        c = cmp(type(v1), type(v2))
    else:
        c = cmp(v1, v2)
    #print 'Comparing %s %s = %s' % (v1, v2, c)
    return c

def calculateBins(variants, field, minValue, maxValue, partitions):
    values = map(lambda x: x.getField(field), variants)
    return calculateBinsForValues(values, field, minValue, maxValue, partitions)

def calculateBinsForValues(values, field, minValue, maxValue, partitions):
    sortedValues = sorted(values)
    captureFieldRangeForPrinting(field, sortedValues)
    
    targetBinSize = len(values) / (1.0*partitions)
    #print sortedValues
    uniqBins = groupby(sortedValues)
    binsAndSizes = map(lambda x: [x[0], len(list(x[1]))], uniqBins)
    #print 'BINS AND SIZES', binsAndSizes

    def bin2Break(bin): return [bin[0], bin[0], bin[1]]
    bins = [bin2Break(binsAndSizes[0])]
    for bin in binsAndSizes[1:]:
        #print '  Breaks', bins
        #print '  current bin', bin
        curLeft = bins[-1][0]
        curSize = bin[1]
        prevSize = bins[-1][2]
        #print curSize, prevSize
        if curSize + prevSize > targetBinSize or (not isNumber(curLeft) and isNumber(bin[0])):
            #print '     => appending', bin2Break(bin)
            bins.append(bin2Break(bin))
        else:
            bins[-1][1] = bin[0]
            bins[-1][2] += curSize

    #print 'Returning ', bins
    #sys.exit(1)
    return bins

def fieldRange(variants, field):
    values = map(lambda v: v.getField(field), variants)
    minValue = min(values)
    maxValue = max(values)
    #rangeValue = maxValue - minValue
    bins = calculateBins(variants, field, minValue, maxValue, OPTIONS.partitions)
    validateBins(bins)
    return minValue, maxValue, bins

def validateBins(bins):
    #print 'Bins are', bins
    for left1, right1, count1 in bins:
        for left2, right2, count2 in bins:
            def contains2(x):
                return left2 < x and x < right2

            if left1 <> left2 and right1 <> right2:
                if None in [left1, left2, right1, right2]:
                    pass # we're ok
                elif contains2(left1) or contains2(right2):
                    raise Exception("Bad bins", left1, right1, left2, right2)

def printFieldQualHeader():
    more = ""
    if TRUTH_CALLS <> None:
        more = CallCmp.HEADER
    def p(stream):
        if stream <> None:
            print >> stream, '  %20s %20s         left        right %15s nVariants  nNovels titv titvNovels  dbSNP  Q' % ("category", "field", "qRange"), more
    p(sys.stdout)
    p(RECAL_LOG)
    
def printFieldQual( category, field, cc, variants ):
    more = ""
    if TRUTH_CALLS <> None:
        callComparison, theseFPs = sensitivitySpecificity(variants, TRUTH_CALLS)
        more = str(callComparison)
    novels = selectVariants(variants, VCFRecord.isNovel)
    def p(stream):
        if stream <> None:
            print >> stream, '  %20s %20s %15s %15s  %8d %8d %2.2f       %2.2f  %3.2f %3d' % (category, field, binString(field, cc.featureRange), cc.qualRange.dashedString(), len(variants), len(novels), titv(variants), titv(novels), dbSNPRate(variants) * 100, phredScale(cc.FPRate)), more
    p(sys.stdout)
    p(RECAL_LOG)

FIELD_RANGES = dict()
def captureFieldRangeForPrinting(field, sortedValues):
    """Finds the minimum float value in sortedValues for convenience printing in recal.log"""
    #print sortedValues
    floatValues = filter(isNumber, sortedValues)
    if floatValues <> []:
        FIELD_RANGES[field] = floatValues[0]
        #print 'Setting field range to', field, FIELD_RANGES[field]


def isNumber(x):
    return isinstance(x, (int, long, float))

def binString(field, cc):
    epsilon = 1e-2
    left, right = cc.left, cc.right
    leftStr = str(left)
    rightStr = "%5s" % str(right)
    if OPTIONS.plottableNones and not isNumber(left) and not isNumber(right) and field in FIELD_RANGES:
        left = right = FIELD_RANGES[field] - epsilon
    if OPTIONS.plottableNones and not isNumber(left) and isNumber(right):
        left = right - epsilon        
    if isNumber(left): leftStr = "%.4f" % left
    if isNumber(right): rightStr = "%.4f" % right
    return '%12s %12s' % (leftStr, rightStr)

#
#
#
def recalibrateCalls(variants, fields, callCovariates):
    def phred(v): return phredScale(v)
    
    for variant in variants:
        recalCall = RecalibratedCall(variant, fields) 
        originalQual = variant.getField('QUAL') 

        for callCovariate in callCovariates:
            if callCovariate.containsVariant(variant):
                FPR = callCovariate.getFPRate()
                recalCall.recalFeature(callCovariate.getFeature(), FPR)
                recalCall.call.setField(callCovariate.getCovariateField(), phred(FPR))

        #recalCall.call.setField('QUAL', phred(recalCall.jointFPErrorRate())) 
        recalCall.call.setField('QUAL', phred(recalCall.jointFPErrorRate())) 
        recalCall.call.setField('OQ', originalQual)
        #print 'recalibrating', variant.getLoc()
        #print '  =>',  variant
        yield recalCall.call
    
#
#
#
def optimizeCalls(variants, fields, titvTarget):
    callCovariates = calibrateFeatures(variants, fields, titvTarget, category = "covariates", useBreaks=True)
    recalCalls = recalibrateCalls(variants, fields, callCovariates)
    return recalCalls, callCovariates

def printCallQuals(field, recalCalls, titvTarget, info = ""):
    print '--------------------------------------------------------------------------------'
    print info
    calibrateFeatures(recalCalls, [field], titvTarget, printCall = True, cumulative = False, forcePrint = True, prefix = "OPT-", printHeader = False, category = "optimized-calls" )
    print 'Cumulative'
    calibrateFeatures(recalCalls, [field], titvTarget, printCall = True, cumulative = True, forcePrint = True, prefix = "OPTCUM-", printHeader = False, category = "optimized-calls" )

def all( p, l ):
    for elt in l:
        if not p(elt): return False
    return True


def mapVariantBins(variants, field, cumulative = False, breakQuals = [Range()]):
    minValue, maxValue, featureBins = fieldRange(variants, field)

    #print 'BREAKQuals', breakQuals[0]
    bins = [(x,y) for x in featureBins for y in breakQuals]
    #print 'BINS', bins

    def variantsInBin(featureBin, qualRange):
        right = featureBin[1]
        if cumulative: 
            right = Range.ANY
        cc = CallCovariate(field, Range(featureBin[0], right), qualRange)

        return cc, selectVariants(variants, lambda v: cc.containsVariant(v))
        
    #sys.exit(1)    
    return starmap( variantsInBin, bins )

def qBreaksRanges(qBreaks, useBreaks):
    if qBreaks == None or not useBreaks:
        return [Range()]        # include everything in a single range
    else:
        breaks = map(float, qBreaks.split(','))
        return map(lambda x, y: Range(x,y, rightOpen = True), chain([Range.ANY], breaks), chain(breaks, [Range.ANY]))

def calibrateFeatures(variants, fields, titvTarget, printCall = False, cumulative = False, forcePrint = False, prefix = '', printHeader = True, category = None, useBreaks = False ):
    covariates = []    

    if printHeader: printFieldQualHeader()
    for field in fields:
        if DEBUG: print 'Optimizing field', field
        
        titv, FPRate = titvFPRateEstimate(variants, titvTarget)
        #print 'Overall FRRate:', FPRate, nErrors, phredScale(FPRate)

        for cc, selectedVariants in mapVariantBins(variants, field, cumulative = cumulative, breakQuals = qBreaksRanges(OPTIONS.QBreaks, useBreaks and field <> 'QUAL')):
            #print 'CC', cc, field, useBreaks
            if len(selectedVariants) > max(OPTIONS.minVariantsPerBin,1) or forcePrint:
                titv, FPRate = titvFPRateEstimate(selectedVariants, titvTarget)
                #dbsnp = dbSNPRate(selectedVariants)
                cc.FPRate = FPRate
                covariates.append(cc)
                printFieldQual( category, prefix + field, cc, selectedVariants )
            else:
                print 'Not calibrating bin', cc, 'because it contains too few variants:', len(selectedVariants)

    return covariates

class CallCmp:
    def __init__(self, nTP, nFP, nFN):
        self.nTP = nTP
        self.nFP = nFP
        self.nFN = nFN
    
#    def FPRate(self):
#        return (1.0*self.nFP) / max(self.nTP + self.nFP, 1)

    def FNRate(self):
        return (1.0*self.nFN) / max(self.nTP + self.nFN, 1)

    def sensitivity(self):
        # = TP / (TP + FN)
        return (1.0*self.nTP) / max( self.nTP + self.nFN,1 )

    def PPV(self):
        # = TP / (TP + FP)
        return (1.0*self.nTP) / max( self.nTP + self.nFP, 1 )
    
    HEADER = "TP    FP    FN  FNRate  Sensitivity PPV"
    
    def __str__(self):
        return '%6d %6d %6d %.3f %.3f %.3f' % (self.nTP, self.nFP, self.nFN, self.FNRate(), self.sensitivity(), self.PPV())

def variantInTruth(variant, truth):
    if variant.getLoc() in truth:
        return truth[variant.getLoc()]
    else:
        return False

def isVariantInSample(t, sample):
    #print "isVariantInSample", t.getLoc(), t.getField(sample), x
    return t.getField(sample) <> "0/0"

def variantsInTruth(truth):
    # fixme
    return len(filter(lambda x: isVariantInSample(x, OPTIONS.useSample), truth))
    
def sensitivitySpecificity(variants, truth):
    nTP, nFP = 0, 0
    FPs = []
    for variant in variants:
        t = variantInTruth(variant, truth)

        isTP, isFP = False, False
        if OPTIONS.useSample or OPTIONS.onlyAtTruth:
            if t: # we have a site
                isTP = (isVariantInSample(t, OPTIONS.useSample) and t.isVariant()) or (not isVariantInSample(t, OPTIONS.useSample) and not t.isVariant())
                isFP = not isTP
        else:
            isTP = t
            isFP = not t

        #if variant.getLoc() == "1:867694":
        #    print variant, 'T: [', t, '] isTP, isFP', isTP, isFP

        if isTP:
            t.setField("FN", 0)
            variant.setField("TP", 1)
            nTP += 1
        elif isFP:
            nFP += 1
            variant.setField("TP", 0)
            #print t, variant, "is a FP!"
            FPs.append(variant)
    nRef = 0 # len(filter(lambda x: not x.isVariant(), truth.itervalues()))
    nFN = variantsInTruth(truth.itervalues()) - nTP - nRef
    #print 'nRef', nTP, nFP, nFN, nRef
    return CallCmp(nTP, nFP, nFN), FPs

def markTruth(calls):
    if not OPTIONS.useSample:
        for variant in calls.itervalues(): 
            variant.setField("TP", 0) # set the TP field to 0

def compareCalls(calls, truthCalls):
    #markTruth(truthCalls)
    
    def compare1(name, cumulative):
        for field in ["QUAL", "OQ"]:
            for cc, selectedVariants in mapVariantBins(calls, field, cumulative = cumulative):
                #print selectedVariants[0]
                printFieldQual("truth-comparison-" + name, field, cc, selectedVariants )
    
    print 'PER BIN nCalls=', len(calls)
    compare1('per-bin', False)

    print 'CUMULATIVE nCalls=', len(calls)
    compare1('cum', True)
    
def randomSplit(l, pLeft):
    def keep(elt, p):
        if p < pLeft:
            return elt, None
        else:
            return None, elt
    data = [keep(elt, p) for elt, p in zip(l, map(lambda x: random.random(), l))]
    def get(i): return filter(lambda x: x <> None, [x[i] for x in data])
    return get(0), get(1)

def setup():
    global OPTIONS, header
    usage = "usage: %prog files.list [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--f", dest="fields",
                        type='string', default="QUAL",
                        help="Comma-separated list of fields (either in the VCF columns of as INFO keys) to use during optimization [default: %default]")
    parser.add_option("-t", "--truth", dest="truth",
                        type='string', default=None,
                        help="VCF formated truth file.  If provided, the script will compare the input calls with the truth calls.  It also emits calls tagged as TP and a separate file of FP calls")
    parser.add_option("-l", "--recalLog", dest="recalLog",
                        type='string', default="recal.log",
                        help="VCF formated truth file.  If provided, the script will compare the input calls with the truth calls.  It also emits calls tagged as TP and a separate file of FP calls")
    parser.add_option("-u", "--unfiltered", dest="unfiltered",
                        action='store_true', default=False,
                        help="If provided, unfiltered calls will be used in comparisons [default: %default]")
    parser.add_option("", "--plottable", dest="plottableNones",
                        action='store_true', default=False,
                        help="If provided, will generate fake plottable points for annotations with None values -- doesn't effect the behavior of the system just makes it easy to plot outputs [default: %default]")
    parser.add_option("", "--onlyAtTruth", dest="onlyAtTruth",
                        action='store_true', default=False,
                        help="If provided, we only consider TP/FP/FN at truth sites[default: %default]")
    parser.add_option("-p", "--partitions", dest="partitions",
                        type='int', default=25,
                        help="Number of partitions to use for each feature.  Don't use so many that the number of variants per bin is very low. [default: %default]")
    parser.add_option("", "--maxRecordsForCovariates", dest="maxRecordsForCovariates",
                        type='int', default=2000000,
                        help="Derive covariate information from up to this many VCF records.  For files with more than this number of records, the system downsamples the reads [default: %default]")
    parser.add_option("-m", "--minVariantsPerBin", dest="minVariantsPerBin",
                       type='int', default=10,
                        help="Don't include any covariates with fewer than this number of variants in the bin, if such a thing happens.  NEEDS TO BE FIXED")
    parser.add_option("-M", "--maxRecords", dest="maxRecords",
                       type='int', default=None,
                        help="Maximum number of input VCF records to process, if provided.  Default is all records")
    parser.add_option("-q", "--qMin", dest="minQScore",
                        type='int', default=-1,
                        help="The minimum Q score of the initial SNP list to consider for selection [default: %default]")
    parser.add_option("-Q", "--qMax", dest="maxQScore",
                        type='int', default=1000000,
                        help="The maximum Q score allowed for both a single covariate and the overall QUAL score [default: %default]")
    parser.add_option("", "--QBreaks", dest="QBreaks",
                        type='string', default=None,
                        help="Breaks in QUAL for generating covarites [default: %default]")
    parser.add_option("", "--maxQScoreForCovariate", dest="maxQScoreForCovariate",
                        type='int', default=60,
                        help="The maximum Q score allowed for both a single covariate and the overall QUAL score [default: %default]")
    parser.add_option("-o", "--outputVCF", dest="outputVCF",
                        type='string', default="recal.vcf",
                        help="If provided, a VCF file will be written out to this file [default: %default]")
    parser.add_option("", "--FNoutputVCF", dest="FNoutputVCF",
                        type='string', default=None,
                        help="If provided, VCF file will be written out to this file [default: %default]")
    parser.add_option("", "--titv", dest="titvTarget",
                        type='float', default=None,
                        help="If provided, we will optimize calls to the targeted ti/tv rather than that calculated from known calls [default: %default]")
    parser.add_option("-b", "--bootstrap", dest="bootStrap",
                       type='float', default=None,
                       help="If provided, the % of the calls used to generate the recalibration tables. [default: %default]")
    parser.add_option("-k", "--skip", dest="skip",
                       type=int, default=None,
                       help="If provided, we'll only take every nth record [default: %default]")
    parser.add_option("-s", "--useSample", dest="useSample",
                        type='string', default=False,
                        help="If provided, we will examine sample genotypes for this sample, and consider TP/FP/FN in the truth conditional on sample genotypes [default: %default]")
    parser.add_option("-r", "--dontRecalibrate", dest="dontRecalibrate",
                        action='store_true', default=False,
                        help="If provided, we will not actually do anything to the calls, they will just be assessed [default: %default]")
   
    (OPTIONS, args) = parser.parse_args()
    if len(args) > 2:
        parser.error("incorrect number of arguments")
    return args

def assessCalls(file):
    print 'Counting records in VCF', file
    numberOfRecords = 1
    downsampleFraction = 1
    if OPTIONS.maxRecords <> None:
        numberOfRecords = quickCountRecords(open(file))
        if OPTIONS.maxRecords < numberOfRecords:
            numberOfRecords = OPTIONS.maxRecords
        downsampleFraction = min(float(OPTIONS.maxRecordsForCovariates) / numberOfRecords, 1)
    #print 'Reading variants', OPTIONS.skip, downsampleFraction
    header, allCalls = readVariants(file, OPTIONS.maxRecords, downsampleFraction=downsampleFraction, minQScore = OPTIONS.minQScore, maxQScore = OPTIONS.maxQScore, skip = OPTIONS.skip)
    allCalls = list(allCalls)
    print 'Number of VCF records', numberOfRecords, ', max number of records for covariates is', OPTIONS.maxRecordsForCovariates, 'so keeping', downsampleFraction * 100, '% of records'
    print 'Number of selected VCF records', len(allCalls)
    
    titvtarget = OPTIONS.titvTarget
    if titvtarget == None:
        titvtarget = titv(selectVariants(allCalls, VCFRecord.isKnown))
    print 'Ti/Tv all  ', titv(allCalls)
    print 'Ti/Tv known', titv(selectVariants(allCalls, VCFRecord.isKnown))
    print 'Ti/Tv novel', titv(selectVariants(allCalls, VCFRecord.isNovel))
    
    return header, allCalls, titvtarget

def determineCovariates(allCalls, titvtarget, fields):
    if OPTIONS.bootStrap:
        callsToOptimize, recalEvalCalls = randomSplit(allCalls, OPTIONS.bootStrap)
    else:
        callsToOptimize = allCalls 

    recalOptCalls, covariates = optimizeCalls(callsToOptimize, fields, titvtarget)
    printCallQuals("QUAL", list(recalOptCalls), titvtarget, 'OPTIMIZED CALLS')

    if OPTIONS.bootStrap:
        recalibatedEvalCalls = recalibrateCalls(recalEvalCalls, fields, covariates)
        printCallQuals("QUAL", list(recalibatedEvalCalls), titvtarget, 'BOOTSTRAP EVAL CALLS')

    return covariates

def writeRecalibratedCalls(file, header, calls):
    if file:
        f = open(file, 'w')
        #print 'HEADER', header
        i = 0
        for line in formatVCF(header, calls):
            if i % 10000 == 0: print 'writing VCF record', i
            i += 1
            print >> f, line
        f.close()

def readTruth(truthVCF):
    print 'Reading truth file', truthVCF
    rawTruth = list(readVariants(truthVCF, maxRecords = None, decodeAll = True, mustBeVariant = True)[1])
    truth = dict( [[v.getLoc(), v] for v in rawTruth])
    print 'Number of raw and passing filter truth calls', len(rawTruth), len(truth)
    return truth

def evaluateTruth(header, callVCF, truth, truthVCF):
    print 'Reading variants back in from', callVCF
    header, calls = readVariants(callVCF)
    calls = list(calls)
    
    print '--------------------------------------------------------------------------------'
    print 'Comparing calls to truth', truthVCF
    print ''

    compareCalls(calls, truth)

    writeRecalibratedCalls(callVCF, header, calls)

    def isFN(v):
        return isVariantInSample(v, OPTIONS.useSample) and not v.hasField("FN")

    if truth <> None and OPTIONS.FNoutputVCF:
        f = open(OPTIONS.FNoutputVCF, 'w')
        #print 'HEADER', header
        for line in formatVCF(header, filter( isFN, truth.itervalues())):
            print >> f, line
        f.close()

TRUTH_CALLS = None
RECAL_LOG = None
def main():
    global TRUTH_CALLS, RECAL_LOG
    
    args = setup()
    fields = OPTIONS.fields.split(',')

    truthVCF = None
    #print("LENGTH OF ARGS "+str(len(args)))

    if OPTIONS.truth <> None:
        truthVCF = OPTIONS.truth
        TRUTH_CALLS = readTruth(truthVCF)

    if OPTIONS.recalLog <> None:
        RECAL_LOG = open(OPTIONS.recalLog, "w") 
        print >> RECAL_LOG, "# optimized vcf", args[0]
        print >> RECAL_LOG, "# truth vcf", truthVCF
        for key, value in OPTIONS.__dict__.iteritems():
            print >> RECAL_LOG, '#', key, value

    header, allCalls, titvTarget = assessCalls(args[0])
    if not OPTIONS.dontRecalibrate: 
        covariates = determineCovariates(allCalls, titvTarget, fields)
        print OPTIONS
        header, callsToRecalibate = readVariants(args[0], OPTIONS.maxRecords, minQScore = OPTIONS.minQScore, maxQScore = OPTIONS.maxQScore, skip = OPTIONS.skip)
        RecalibratedCalls = recalibrateCalls(callsToRecalibate, fields, covariates)
        writeRecalibratedCalls(OPTIONS.outputVCF, header, RecalibratedCalls)
    else:
        printFieldQualHeader()
        printCallQuals("QUAL", allCalls, titvTarget)
        OPTIONS.outputVCF = args[0]

    if truthVCF <> None:
        evaluateTruth(header, OPTIONS.outputVCF, TRUTH_CALLS, truthVCF)


PROFILE = False
if __name__ == "__main__":
    if PROFILE:
        import cProfile
        cProfile.run('main()', 'fooprof')
        import pstats
        p = pstats.Stats('fooprof')
        p.sort_stats('cumulative').print_stats(10)
        p.sort_stats('time').print_stats(10)
        p.sort_stats('time', 'cum').print_stats(.5, 'init')
    else:
        main()
