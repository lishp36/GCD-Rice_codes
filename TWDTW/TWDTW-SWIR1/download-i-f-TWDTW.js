// download-i-f.js
// 合成影像，插值、滤波后并下载
// i stand for interpolation
// f  stand for filter
//

/**
 * 对去云后的影像做线性插值
 * @param {ee.ImageCollection} imgCollection 输入的影像集
 * @param range 插值使用的影像（前后多少幅影像）
 * @param nBands 用到的波段数
 * @param bands 输入的影像集里单幅影像的波段名。例：["swir1", "swir2"]
 * @returns 插值后的影像
 */
exports.linear_interpolation = function (imgCollection, range, nBands, bands) {
    var size = imgCollection.size();
    var imgList = imgCollection.toList(size);
    var interpolated = ee.ImageCollection(ee.List.sequence(
      range - 1, nBands - range, 1
    ).map(function (i) {
      var i = ee.Number(i)
      var before = ee.ImageCollection.fromImages(imgList.slice(i.subtract(range - 1), i))
        .filter(ee.Filter.gt('length', 0)).mosaic();
      var after = ee.ImageCollection.fromImages(imgList.slice(i.add(1), i.add(range)).reverse())
        .filter(ee.Filter.gt('length', 0)).mosaic();
      var boforeY = before.select(bands);
      var beforedoy = before.select("index");
      var afterY = after.select(bands);
      var afterdoy = after.select("index");
      var targetImg = ee.Image(imgList.get(i));
      var currentdoy = ee.Image.constant(targetImg.get("index")).float();
      var Y = afterY.subtract(boforeY).divide(afterdoy.subtract(beforedoy))
        .multiply(currentdoy.subtract(beforedoy)).add(boforeY);
      var filledImage = ee.Image(ee.Algorithms.If({
        condition: ee.Number(targetImg.get('length')).gt(0),
        trueCase: targetImg.select(bands).unmask(Y),
        falseCase: Y
      }));
      return filledImage.unmask(0).clip(studyPlace)
        .set('system:time_start', targetImg.get('system:time_start'), "index", targetImg.get("index")); // can not simply copy all properties of imgCollection
    }));
    return interpolated;
  };
  
  /**
   * 对影像做SG滤波
   * @param {ee.ImageCollection} collection 输入的影像集
   * @param nBands 输入的影像集里单幅影像的所需的波段数
   * @param iBand 需要插值的波段的序号
   * @param bandName 需要插值的波段名（一次只能对某一波段滤波）。例："swir1"
   * @returns 滤波后的影像
   */
  exports.SGFilter = function (collection, nBands, iBand, bandName) {
    // Step 2: Set up Savitzky-Golay smoothing
    var window_size = 5;
    var half_window = (window_size - 1) / 2;
  
    // Define the axes of variation in the collection array.
    var imageAxis = 0;
    var bandAxis = 1;
  
    // Set polynomial order
    // var order = 3;
    // var coeffFlattener = [["constant", "x", "x2", "x3"]];
    // var indepSelectors = ["constant", "t", "t2", "t3"];
  
    // Change to order = 2 as follows:
    var order = 2;
    var coeffFlattener = [["constant", "x", "x2"]];
    var indepSelectors = ["constant", "t", "t2"];
  
    // Convert the collection to an array.
    var array = collection.toArray();
  
    // Solve
    function getLocalFit(i) {
      // Get a slice corresponding to the window_size of the SG smoother
      var subarray = array.arraySlice(imageAxis, ee.Number(i).int(), ee.Number(i).add(window_size).int());
      var predictors = subarray.arraySlice(bandAxis, nBands, nBands + order + 1);
      var response = subarray.arraySlice(bandAxis, iBand, iBand + 1);
      var coeff = predictors.matrixSolve(response);
  
      coeff = coeff.arrayProject([0]).arrayFlatten(coeffFlattener);
      return coeff;
    }
  
    // For the remainder, use collection as a list of images
    var collectionList = collection.toList(collection.size());
    var runLength = ee.List.sequence(0, collectionList.size().subtract(window_size));
  
    // Run the SG solver over the series, and return the smoothed image version
    var sg_series = runLength.map(function (i) {
      var ref = ee.Image(collectionList.get(ee.Number(i).add(half_window)));
      return getLocalFit(i)
        .multiply(ref.select(indepSelectors)).reduce(ee.Reducer.sum())
        .rename(bandName).copyProperties(ref)
        ;
    });
    return sg_series;
  };
  
  /**
   * 指定时间段和波段，获取 SAR 和 光学（S2, L7, L8）影像
   * @param dStart1 SAR 开始日期
   * @param dEnd1 SAR 结束日期
   * @param dStart2 光学开始日期
   * @param dEnd2 光学结束日期
   * @param {ee.Feature} studyPlace 研究区域的范围
   * @param {Function} pre 预处理函数
   * @param dStep1 SAR 合成（日/月数）
   * @param unit1 SAR 合成单位 day/month
   * @param dStep2 光学合成（日/月数）
   * @param unit2 光学合成单位 day/month
   * @param {ee.Image} mask 耕地掩膜
   * @param {List} band1 SAR 波段
   * @param {List} band2 光学波段
   * @returns 合成的影像（时间合成到波段）
   */
  exports.getCollection = function (
    yr, dStart1, dEnd1, dStart2, dEnd2,
    studyPlace, pre, dStep1, unit1, dStep2, unit2, mask, band1, band2, description
  ) {
    var range = 5;
    var dateFilter2 = ee.Filter.date(
      dStart2.advance(-(range + 2) * dStep2, unit2),
      dEnd2.advance((range + 2) * dStep2, unit2)
    );
    var S2 = ee.ImageCollection('COPERNICUS/S2')
      .filter(dateFilter2).filterBounds(studyPlace);
    var S2C = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
      .filter(dateFilter2).filterBounds(studyPlace);
    var collection_S2 = pre.indexJoin(S2, S2C, 'cloud_probability')
      .map(pre.scaleFactorsS2).map(pre.maskImage_S2)//.map(pre.rmCloudS2)
      .map(pre.renameBandsS2)
      .map(pre.getNDVI).map(pre.getLSWI).map(pre.getREP).map(pre.getGCVI)
      ;
    var collection_L5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
      .merge(ee.ImageCollection("LANDSAT/LT05/C02/T2_L2"))
      .filter(dateFilter2).filterBounds(studyPlace)
      .map(pre.scaleFactorsLC2).map(pre.renameBandsL7C2).map(pre.rmCloudLC2)
      .map(pre.getNDVI).map(pre.getGCVI)
      ;
    var collection_L7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
      .merge(ee.ImageCollection("LANDSAT/LE07/C02/T2_L2"))
      .filter(dateFilter2).filterBounds(studyPlace)
      .map(pre.scaleFactorsLC2).map(pre.renameBandsL7C2).map(pre.rmCloudLC2)
      .map(pre.getNDVI).map(pre.getGCVI)
      ;
    var collection_L8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
      .merge(ee.ImageCollection("LANDSAT/LC08/C02/T2_L2"))
      .filter(dateFilter2).filterBounds(studyPlace)
      .map(pre.scaleFactorsLC2).map(pre.renameBandsL8C2).map(pre.rmCloudLC2)
      .map(pre.getNDVI).map(pre.getGCVI)
      ;
    var collection_L9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
      .merge(ee.ImageCollection("LANDSAT/LC09/C02/T2_L2"))
      .filter(dateFilter2).filterBounds(studyPlace)
      .map(pre.scaleFactorsLC2).map(pre.renameBandsL8C2).map(pre.rmCloudLC2)
      .map(pre.getNDVI).map(pre.getGCVI)
      ;
    var collection_opt =
      // collection_S2
      collection_L8
      .merge(collection_L5)
      // .merge(collection_L7) // 1999 - 2001, 2012除了这几年，其余年份注释， 线性扫描器坏
      // .merge(collection_L8)
      .merge(collection_L9)
      .select(band2)
      ;
  
    var dateFilter1 = ee.Filter.date(
      dStart1.advance(-(range + 2) * dStep1, unit1),
      dEnd1.advance((range + 2) * dStep1, unit1)
    );
    var collection_S1 = ee.ImageCollection("COPERNICUS/S1_GRD")
      .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
      .filter(ee.Filter.eq('instrumentMode', 'IW'))
      // .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
      .filter(dateFilter1).filterBounds(studyPlace)
      ;
  
    // composte
    var size2 = dEnd2.difference(dStart2, unit2).divide(dStep2).toInt().getInfo();
    var list2 = ee.List.sequence(-range - 1, size2 + range, 1);
    var collection2 = ee.ImageCollection(list2.map(function (d) {
      var d_start = dStart2.advance(ee.Number(d).multiply(dStep2), unit2);
      var d_end = d_start.advance(dStep2, unit2);
      var collection_opt_ = collection_opt.filterDate(d_start, d_end)
        .select(band2).median()
        ;
      var bandLength2 = collection_opt_.bandNames().length();
      var mask2 = ee.Algorithms.If({
        condition: ee.Number(bandLength2).gt(0),
        trueCase: collection_opt_.select(0).mask(),
        falseCase: ee.Image(0).clip(studyPlace)
      });
      return collection_opt_
        .addBands(ee.Image.constant(d).rename("index").float())
        .updateMask(mask2)
        .set('system:time_start', d_start.millis())
        
        .set("index", d)
        .set('length', bandLength2)
        .addBands(ee.Image(mask2).rename("mask"))
        ;
    }));
  
    var opt_valid0 = collection2.select("mask")
      .filter(ee.Filter.gte("index", 0))
      .filter(ee.Filter.lte("index", size2-1));
    var opt_valid = opt_valid0.toBands().clip(studyPlace).updateMask(mask).toInt8();
    var opt_num = opt_valid0.sum().clip(studyPlace).updateMask(mask).toInt16();
    // Map.addLayer(opt_num, {min: 0, max: size2-1}, "S2_num");
    // Map.addLayer(collection2.select(band2[0])
    //   .filter(ee.Filter.gte("index", 0))
    //   .filter(ee.Filter.lte("index", size2-1)).toBands(),
    //   {}, band2[0]
    // );
    
    Export.image.toDrive({
      image: opt_valid,
      description: description + "-valid",
      scale: 30,
      fileNamePrefix: description + "-valid-WGS84",
      maxPixels: 1e13,
      region: studyPlace,
      folder: "",
      crs: "EPSG:4326",
    });

    var size1 = dEnd1.difference(dStart1, unit1).divide(dStep1).toInt().getInfo();
    var list1 = ee.List.sequence(-range - 1, size1 + range, 1);
    var collectionS1 = ee.ImageCollection(list1.map(function (d) {
      var d_start = dStart1.advance(ee.Number(d).multiply(dStep1), unit1);
      var d_end = d_start.advance(dStep1, unit1);
      var collection_S1_ = collection_S1.filterDate(d_start, d_end)
        .median().select(band1)
        ;
      var bandLength1 = collection_S1_.bandNames().length();
      var mask1 = ee.Algorithms.If({
        condition: ee.Number(bandLength1).gt(0),
        trueCase: collection_S1_.select(0).mask(),
        falseCase: ee.Image(0).clip(studyPlace).selfMask()
      });
      return collection_S1_
        .addBands(ee.Image.constant(d).rename("index").float())
        .updateMask(mask1)
        .set('system:time_start', d_start.millis())
        .set("index", d)
        .set('length', bandLength1)
        ;
    }));
  
    var collection2_ = exports.linear_interpolation(
      collection2, range, size2 + 2 * range + 2, band2
    );
    collection2_ = collection2_.map(function (img) {
      var t = ee.Image.constant(img.get("index")).float();
      return img
        .addBands(ee.Image(1).toFloat().rename("constant"))
        .addBands(t.rename("t"))
        .addBands(t.pow(ee.Image(2)).rename("t2"))
        // .addBands(t.pow(ee.Image(3)).rename("t3"))
        ;
    });
    var collectionS1_ = exports.linear_interpolation(
      collectionS1, range, size1 + 2 * range + 2, band1
    );
    collectionS1_ = collectionS1_.map(function (img) {
      var t = ee.Image.constant(img.get("index")).float();
      return img
        .addBands(ee.Image(1).toFloat().rename("constant"))
        .addBands(t.rename("t"))
        .addBands(t.pow(ee.Image(2)).rename("t2"))
        // .addBands(t.pow(ee.Image(3)).rename("t3"))
        ;
    });
  
    if (band1.length > 0) {
      var result = ee.ImageCollection(exports.SGFilter(
        collectionS1_, band1.length, 0, band1[0]
      )).toBands();
      for (var j = 1; j < band1.length; j++) {
        result = result.addBands(ee.ImageCollection(exports.SGFilter(
          collectionS1_, band1.length, j, band1[j]
        )).toBands());
      }
      for (var j = 0; j < band2.length; j++) {
        result = result.addBands(ee.ImageCollection(exports.SGFilter(
          collection2_, band2.length, j, band2[j]
        )).toBands());
      }
      return result.clip(studyPlace).updateMask(mask);
    } else {
      var result = ee.ImageCollection(exports.SGFilter(
        collection2_, band2.length, 0, band2[0]
      )).toBands();
      for (var j = 1; j < band2.length; j++) {
        result = result.addBands(ee.ImageCollection(exports.SGFilter(
          collection2_, band2.length, j, band2[j]
        )).toBands());
      }
      for (var j = 0; j < band1.length; j++) {
        result = result.addBands(ee.ImageCollection(exports.SGFilter(
          collectionS1_, band1.length, j, band1[j]
        )).toBands());
      }
      return result.clip(studyPlace).updateMask(mask);
    }
  
  };
  
  var pre = require("");
  var At = ""; var province = At;
  var chinaProvinces = ee.FeatureCollection("");
  var studyPlace = chinaProvinces.filterMetadata("", "equals", "");

  
  var band1 = [
    // "VV",
    // "VH",
  ];
  var band2 = [
    // "blue",
    // "green",
    // "red",
    // "re2",
    // "nir",
    // "REP",
    "swir1",
    // "swir2",
    // "NDVI",
    // "LSWI",
    // "GCVI",
  ];
  Map.addLayer(studyPlace, {}, "study place");
  
  for (var yr = 2016; yr <= 2022; yr = yr + 1) {
    var mask0 = ee.Image(
      "users/shinzoqchiuq/fromglc10_2017v01/" + province + "-cropland-WGS84"
    );

    var mask = ee.Image.constant(1);
    // Map.addLayer(mask.selfMask());
  

    var start1 = ""; var end1 = ""; var dStep1 = 12; var unit1 = "day";
    var dStart1 = ee.Date.parse("yyyyDDD", yr + start1);
    var dEnd1 = ee.Date.parse("yyyyDDD", yr + end1)//.advance(dStep1, unit1);
  
    var start2 = ""; var end2 = ""; var dStep2 = 8; var unit2 = "day";
    var dStart2 = ee.Date.parse("yyyyDDD", yr + start2);
    var dEnd2 = ee.Date.parse("yyyyDDD", yr + end2).advance(dStep2, unit2);
  

  
    var size1 = dEnd1.difference(dStart1, unit1).divide(dStep1).toInt().getInfo();
    var n = 0;
    var bands1 = [];
    for (var j = 0; j < band1.length; j++) {
      for (var i = 0; i < size1; i++) {
        bands1[n] = band1[j] + "_" + i.toString(); n += 1;
      }
    }
    var size2 = dEnd2.difference(dStart2, unit2).divide(dStep2).toInt().getInfo();
    var n = 0;
    var bands2 = [];
    for (var j = 0; j < band2.length; j++) {
      for (var i = 0; i < size2; i++) {
        bands2[n] = band2[j] + "_" + i.toString(); n += 1;
      }
    }
    var bands = bands1.concat(bands2);
    print(bands);
    var description = At + "-" + yr;
    if (band1.length > 0) {
      description = description + "-S1";
      for (var j = 0; j < band1.length; j++) {
        description = description + "_" + band1[j];
      }
      description = description + "-" + start1 + "_" + dStep1 + "_" + end1 + "_" + unit1;
    }
    if (band2.length > 0) {
      description = description + "-L5_L7";
      for (var j = 0; j < band2.length; j++) {
        description = description + "_" + band2[j];
      }
      if (dStep2 == "0.5") {
      description = description + "-" + start2 + "_half_" + end2 + "_" + unit2;
      } else {
        description = description + "-" + start2 + "_" + dStep2 + "_" + end2 + "_" + unit2;
      }
    }
    description = description + "-i-f";
  
    var compositedImageClip = exports.getCollection(
      yr, dStart1, dEnd1, dStart2, dEnd2,
      studyPlace, pre, dStep1, unit1, dStep2, unit2, mask, band1, band2, description
    ).rename(bands).subtract(0).multiply(1000).toInt16();
  
    print("compositedImageClip", compositedImageClip);
    Map.addLayer(compositedImageClip, {}, "Image");
  
    Export.image.toDrive({
      image: compositedImageClip,
      description: description,
      scale: 10,
      fileNamePrefix: description + "-WGS84",
      maxPixels: 1e13,
      region: studyPlace,
      folder: "",
      crs: "EPSG:4326",
    });

  }
  