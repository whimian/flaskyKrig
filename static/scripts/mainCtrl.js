angular
    .module('ngYU')
    .controller('mainCtrl', ['$scope', '$sce', 'mainFac', mainCtrl]);

function mainCtrl($scope, $sce, mainFac) {
    $scope.img_src;

    $scope.getResult = function() {
        mainFac.getResult().success(function(data) {
            $scope.img_src = data;
            console.log($scope.img_src);
        }).error(function() {
            console.log(error);
        });
    };

    $scope.trustDangerousSnippet = function() {
        return $sce.trustAsHtml($scope.img_src);
    };
}
