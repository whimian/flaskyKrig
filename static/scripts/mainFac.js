angular
    .module('ngYU')
    .factory('mainFac',  ['$http', mainFac]);

function mainFac($http) {

    function getResult() {
        return $http.get('/cal');
    }

    return {
        getResult: getResult
    }
}
