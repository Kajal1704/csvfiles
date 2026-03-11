// // import React from "react";

// // const ModelComparison = () => {
// //   return (
// //     <div className="min-h-screen bg-gray-50 p-10">

// //       <h1 className="text-3xl font-bold text-gray-800 mb-8">
// //         Model Performance Comparison
// //       </h1>

// //       {/* Comparison Cards */}
// //       <div className="grid md:grid-cols-2 gap-8">

// //         {/* Random Forest */}
// //         <div className="bg-white rounded-2xl shadow-md p-6">
// //           <h2 className="text-xl font-semibold text-teal-700 mb-4">
// //             Random Forest
// //           </h2>

// //           <ul className="text-gray-600 space-y-2">
// //             <li>• Baseline classification model</li>
// //             <li>• Good overall accuracy</li>
// //             <li>• Risk Score identified as most influential feature</li>
// //             <li>• Moderate generalization</li>
// //           </ul>

// //           <div className="h-48 bg-gray-100 rounded-lg mt-6 flex items-center justify-center text-gray-400">
// //             Random Forest Confusion Matrix
// //           </div>
// //         </div>

// //         {/* XGBoost */}
// //         <div className="bg-white rounded-2xl shadow-md p-6">
// //           <h2 className="text-xl font-semibold text-teal-700 mb-4">
// //             XGBoost (Selected Model)
// //           </h2>

// //           <ul className="text-gray-600 space-y-2">
// //             <li>• Higher accuracy (~91%)</li>
// //             <li>• Better generalization</li>
// //             <li>• Strong feature importance ranking</li>
// //             <li>• Selected as final prediction model</li>
// //           </ul>

// //           <div className="h-48 bg-gray-100 rounded-lg mt-6 flex items-center justify-center text-gray-400">
// //             XGBoost Confusion Matrix
// //           </div>
// //         </div>

// //       </div>

// //       {/* Conclusion Section */}
// //       <div className="bg-white rounded-2xl shadow-md p-6 mt-10">
// //         <h2 className="text-xl font-semibold text-gray-800 mb-4">
// //           Final Model Selection
// //         </h2>
// //         <p className="text-gray-600 leading-relaxed">
// //           Based on comparative evaluation, XGBoost demonstrated better
// //           classification accuracy and generalization performance across
// //           disease classes. Therefore, it was selected as the final model
// //           for deployment and SHAP-based explainability.
// //         </p>
// //       </div>

// //     </div>
// //   );
// // };

// // export default ModelComparison;


// import React, { useEffect, useState } from "react";

// const ModelComparison = () => {
//   const [data, setData] = useState(null);

//   useEffect(() => {
//     fetch("http://127.0.0.1:8000/model-comparison")
//       .then(res => res.json())
//       .then(result => setData(result))
//       .catch(err => console.error(err));
//   }, []);

//   if (!data) {
//     return <div className="p-10">Loading model comparison...</div>;
//   }

//   return (
//     <div className="min-h-screen bg-gray-50 px-10 py-10">
//       <h1 className="text-3xl font-bold mb-10">
//         Model Performance Comparison
//       </h1>

//       <div className="grid grid-cols-1 md:grid-cols-2 gap-8">

//         {/* Random Forest */}
//         <div className="bg-white rounded-2xl shadow-md p-6">
//           <h2 className="text-xl font-semibold mb-2 text-green-700">
//             Random Forest
//           </h2>

//           <p className="mb-4 font-medium">
//             Accuracy: {data.rf_accuracy}%
//           </p>

//           <img
//             src={`data:image/png;base64,${data.rf_confusion_matrix}`}
//             alt="RF Confusion Matrix"
//           />
//         </div>

//         {/* XGBoost */}
//         <div className="bg-white rounded-2xl shadow-md p-6">
//           <h2 className="text-xl font-semibold mb-2 text-blue-700">
//             XGBoost (Selected Model)
//           </h2>

//           <p className="mb-4 font-medium">
//             Accuracy: {data.xgb_accuracy}%
//           </p>

//           <img
//             src={`data:image/png;base64,${data.xgb_confusion_matrix}`}
//             alt="XGB Confusion Matrix"
//           />
//         </div>

//       </div>
//     </div>
//   );
// };

// export default ModelComparison;


import { useEffect, useState } from "react";

export default function ModelComparison() {

  const [data, setData] = useState(null);

  useEffect(() => {

    fetch("http://127.0.0.1:8000/model-comparison")
      .then(res => res.json())
      .then(result => {
        setData(result);
      })
      .catch(err => console.error(err));

  }, []);

  if (!data) {
    return (
      <div className="p-10 text-gray-500">
        Loading model comparison...
      </div>
    );
  }

  return (

    <div className="min-h-screen bg-gray-50 px-10 py-10">

      <h1 className="text-3xl font-bold mb-10">
        Model Performance Comparison
      </h1>

      <div className="grid grid-cols-2 gap-8">

        {/* Random Forest */}
        <div className="bg-white rounded-xl shadow-md p-6">

          <h2 className="text-green-700 font-semibold text-lg mb-3">
            Random Forest
          </h2>

          {/* <p className="mb-4">
            Accuracy: {data.rf_accuracy}%
          </p> */}

          <img
            src={`data:image/png;base64,${data.rf_confusion_matrix}`}
            alt="RF Confusion Matrix"
          />

        </div>


        {/* XGBoost */}
        <div className="bg-white rounded-xl shadow-md p-6">

          <h2 className="text-blue-700 font-semibold text-lg mb-3">
            XGBoost (Selected Model)
          </h2>

          {/* <p className="mb-4">
            Accuracy: {data.xgb_accuracy}%
          </p> */}

          <img
            src={`data:image/png;base64,${data.xgb_confusion_matrix}`}
            alt="XGB Confusion Matrix"
          />

        </div>

      </div>

    </div>

  );
}