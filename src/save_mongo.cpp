#include "save_mongo.h"

#include <string>
#include <bsoncxx/json.hpp>
#include <mongocxx/client.hpp>
#include <mongocxx/stdx.hpp>
#include <mongocxx/uri.hpp>
#include <mongocxx/instance.hpp>
#include <bsoncxx/builder/basic/helpers.hpp>
#include <bsoncxx/builder/basic/document.hpp>
#include <bsoncxx/builder/basic/array.hpp>
#include <bsoncxx/document/value.hpp>
#include <bsoncxx/document/view_or_value.hpp>

using namespace std;
using bsoncxx::builder::basic::document;
using bsoncxx::builder::basic::array;
using bsoncxx::builder::basic::kvp;
using bsoncxx::document::view_or_value;
using bsoncxx::document::value;

static mongocxx::instance instance{};

void save_mongo(const size_t iter, const arma::mat& gamma, const arma::vec& beta, const arma::mat& mu, const arma::mat& D, const double rss, const double diff, const double loglik)
{
    mongocxx::uri uri("mongodb://admin:HuYGChenzxc0559@localhost:27017");
    mongocxx::client client(uri);
    mongocxx::database db = client["hgwr"];
    mongocxx::collection coll_coef = db["coef"];
    auto builder = document{};
    builder.append(kvp("iter", (int)iter));
    arma::vec intercept = gamma.col(0) + mu.col(0) + beta(0);
    /// Save intercept
    bsoncxx::builder::basic::array intercept_array_value;
    for (auto &&i : intercept)
    {
        intercept_array_value.append(i);
    }
    builder.append(kvp("intercept", intercept_array_value));
    /// Save gamma
    for (size_t c = 1; c < gamma.n_cols; c++)
    {
        bsoncxx::builder::basic::array array_value;
        for (auto &&i : gamma.col(c))
        {
            array_value.append(i);
        }
        string key = "gamma" + to_string(c);
        builder.append(kvp(key, array_value));
    }
    /// Save mu
    for (size_t c = 1; c < mu.n_cols; c++)
    {
        bsoncxx::builder::basic::array array_value;
        for (auto &&i : mu.col(c))
        {
            array_value.append(i);
        }
        string key = "mu" + to_string(c);
        builder.append(kvp(key, array_value));
    }
    /// Save D
    arma::vec D_trimatl = D(arma::trimatl_ind(arma::size(D)));
    bsoncxx::builder::basic::array Dtrimatl_array_value;
    for (auto &&i : D_trimatl)
    {
        Dtrimatl_array_value.append(i);
    }
    builder.append(kvp("Dtrimatl", Dtrimatl_array_value));
    /// Save beta
    for (size_t i = 1; i < beta.n_elem; i++)
    {
        string key = "beta" + to_string(i);
        builder.append(kvp(key, beta(i)));
    }
    /// Save statistics
    builder.append(kvp("RSS", rss));
    builder.append(kvp("dRSS", diff));
    builder.append(kvp("loglik", loglik));
    bsoncxx::document::value doc = builder.extract();
    /// Insert or update document
    auto filter_builder = document{};
    filter_builder.append(kvp("iter", (int)iter));
    auto filter = filter_builder.extract();
    bsoncxx::stdx::optional<bsoncxx::document::value> maybe_result = coll_coef.find_one(view_or_value(filter));
    if (maybe_result)
    {
        auto update_doc_builder = document{};
        update_doc_builder.append(kvp("$set", doc));
        coll_coef.update_one(view_or_value(filter), view_or_value(update_doc_builder));
    }
    else
        bsoncxx::stdx::optional<mongocxx::result::insert_one> result = coll_coef.insert_one(view_or_value(doc));
}
